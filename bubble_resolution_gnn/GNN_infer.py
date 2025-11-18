import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader

# --------------------------
# Utilities / tokenization
# --------------------------

NUC2ID = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 5}
ID2NUC = {v: k for k, v in NUC2ID.items()}

def seq_to_tokens(seq: str) -> torch.LongTensor:
    return torch.tensor([NUC2ID.get(c, 5) for c in seq], dtype=torch.long)

def _safe_get(d: Dict[str, Any], key: str, default):
    v = d.get(key, default)
    return default if v is None else v

# --------------------------
# Dataset for BubbleGun JSONL
# --------------------------

class BubbleGunBubbleDataset(Dataset):
    """
    Dataset built directly from:
      - BubbleGun JSON/JSONL (with numeric node IDs)
      - Graph GFA (from Rust assembler) for seq + coverage + adjacency

    For each bubble:
      - nodes = inside -a endpoints
      - node features: seq_tokens, KC coverage
      - edges = all graph edges u->v where u,v are in that node set
      - edge features: [len_nodes=1, len_bp=k, cov_min=EC, cov_mean=EC, on_ref=0]
    """
    def __init__(self, bubbles_path: str, gfa_path: str):
        super().__init__()
        self.bubbles_path = bubbles_path
        self.gfa_path = gfa_path

        # --- parse GFA once ---
        self.node_seq: Dict[str, str] = {}
        self.node_cov: Dict[str, float] = {}
        self.adj: Dict[str, List[Tuple[str, float]]] = {}
        self.k: int = 0

        with open(gfa_path, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("H"):
                    continue
                cols = line.split("\t")
                tag0 = cols[0]
                if tag0 == "S" and len(cols) >= 3:
                    nid = cols[1]
                    seq = cols[2]
                    kc = 0
                    for field in cols[3:]:
                        lf = field.lower()
                        if lf.startswith("kc:i:"):
                            try:
                                kc = int(field.split(":")[-1])
                            except Exception:
                                kc = 0
                    self.node_seq[nid] = seq
                    self.node_cov[nid] = float(kc)
                    if self.k == 0:
                        self.k = len(seq)
                elif tag0 == "L" and len(cols) >= 7:
                    u = cols[1]
                    v = cols[3]
                    ec = 0
                    for field in cols[6:]:
                        lf = field.lower()
                        if lf.startswith("ec:i:"):
                            try:
                                ec = int(field.split(":")[-1])
                            except Exception:
                                ec = 0
                    self.adj.setdefault(u, []).append((v, float(ec)))

        if self.k == 0 and self.node_seq:
            # Fallback: take length of any seq
            self.k = len(next(iter(self.node_seq.values())))

        # --- parse BubbleGun JSON/JSONL into bubbles ---
        self.bubbles: List[Dict[str, Any]] = []
        with open(bubbles_path, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    obj = json.loads(line)
                except json.JSONDecodeError:
                    continue

                chains: List[Dict[str, Any]] = []
                if isinstance(obj, dict) and "bubbles" in obj:
                    chains = [obj]
                elif isinstance(obj, dict):
                    chains = [v for v in obj.values() if isinstance(v, dict)]
                elif isinstance(obj, list):
                    chains = [c for c in obj if isinstance(c, dict)]
                else:
                    continue

                for ch in chains:
                    bs = ch.get("bubbles", [])
                    for b in bs:
                        if not isinstance(b, dict):
                            continue
                        ends = b.get("ends", [])
                        inside = b.get("inside", [])
                        if not isinstance(ends, list) or len(ends) != 2:
                            continue
                        start_id = str(ends[0])
                        end_id = str(ends[1])
                        inside_ids = [str(x) for x in inside]
                        self.bubbles.append(
                            {
                                "bubble_id": b.get("id", len(self.bubbles)),
                                "start": start_id,
                                "end": end_id,
                                "inside": inside_ids,
                            }
                        )

    def __len__(self) -> int:
        return len(self.bubbles)

    def __getitem__(self, idx: int) -> Data:
        rec = self.bubbles[idx]
        node_ids = set(rec["inside"])
        node_ids.add(rec["start"])
        node_ids.add(rec["end"])

        # Map node IDs -> local indices and gather seq + cov
        id_to_local: Dict[str, int] = {}
        node_seqs: List[str] = []
        node_covs: List[float] = []

        for nid in node_ids:
            seq = self.node_seq.get(nid)
            if seq is None:
                # if BubbleGun refers to node missing from GFA, skip this node
                continue
            local_idx = len(node_seqs)
            id_to_local[nid] = local_idx
            node_seqs.append(seq)
            cov = self.node_cov.get(nid, 0.0)
            node_covs.append(float(cov))

        n_nodes = len(node_seqs)
        if n_nodes == 0:
            # empty graph, but keep record shape consistent
            data = Data(
                seq_tokens=torch.zeros((0, 1), dtype=torch.long),
                x_cov=torch.zeros((0, 1), dtype=torch.float32),
                edge_index=torch.zeros((2, 0), dtype=torch.long),
                edge_attr=torch.zeros((0, 5), dtype=torch.float32),
                num_nodes=0,
            )
            data.node_seqs = []
            data.bubble_id = rec["bubble_id"]
            return data

        # Build edges restricted to this bubble's node set
        edge_src: List[int] = []
        edge_dst: List[int] = []
        edge_attr: List[List[float]] = []

        for u in node_ids:
            if u not in id_to_local:
                continue
            for (v, ec) in self.adj.get(u, []):
                if v not in id_to_local:
                    continue
                u_idx = id_to_local[u]
                v_idx = id_to_local[v]
                # len_nodes=1, len_bp=k, cov_min=cov_mean=EC, on_ref=0
                edge_src.append(u_idx)
                edge_dst.append(v_idx)
                edge_attr.append([1.0, float(self.k), ec, ec, 0.0])

        if edge_src:
            edge_index = torch.tensor([edge_src, edge_dst], dtype=torch.long)
            edge_attr_t = torch.tensor(edge_attr, dtype=torch.float32)
        else:
            edge_index = torch.zeros((2, 0), dtype=torch.long)
            edge_attr_t = torch.zeros((0, 5), dtype=torch.float32)

        seq_tokens = torch.stack([seq_to_tokens(s) for s in node_seqs], dim=0)
        x_cov = torch.tensor(node_covs, dtype=torch.float32).unsqueeze(1)

        data = Data(
            seq_tokens=seq_tokens,
            x_cov=x_cov,
            edge_index=edge_index,
            edge_attr=edge_attr_t,
            num_nodes=seq_tokens.size(0),
        )
        data.node_seqs = node_seqs
        data.bubble_id = rec["bubble_id"]
        data.k = self.k
        return data


# --------------------------
# Model
# --------------------------

class DenseGCNLayer(nn.Module):
    def __init__(self, in_dim, out_dim):
        super().__init__()
        self.lin = nn.Linear(in_dim, out_dim)

    def forward(self, H, edge_index_local, n_nodes: int):
        if n_nodes == 0:
            return H
        A = H.new_zeros((n_nodes, n_nodes))
        if edge_index_local.numel() > 0:
            src, dst = edge_index_local
            one = torch.ones_like(src, dtype=H.dtype)
            A.index_put_((src, dst), one, accumulate=True)
        idx = torch.arange(n_nodes, device=H.device)
        try:
            A[idx, idx] += 1
        except Exception:
            flat = A.view(-1)
            diag_idx = torch.arange(0, n_nodes*n_nodes, step=n_nodes+1, device=H.device)
            flat.index_add_(0, diag_idx, torch.ones(n_nodes, device=H.device, dtype=H.dtype))
            A = flat.view(n_nodes, n_nodes)
        deg = A.sum(dim=1) + 1e-8
        D_inv_sqrt = deg.pow(-0.5)
        A_norm = (D_inv_sqrt[:, None] * A) * D_inv_sqrt[None, :]
        return A_norm @ self.lin(H)

class HyperbubbleGNN(nn.Module):
    def __init__(
        self,
        vocab_size=6,
        embed_dim=16,
        gcn_hidden=16,
        edge_hidden=16,
        edge_feat_dim=5,
        use_lstm=False,
    ):
        super().__init__()
        self.embed = nn.Embedding(vocab_size, embed_dim, padding_idx=0)
        self.use_lstm = use_lstm
        if use_lstm:
            self.lstm = nn.LSTM(embed_dim, gcn_hidden // 2, batch_first=True, bidirectional=False)
            self.node_in = gcn_hidden + 1
        else:
            self.node_in = embed_dim + 1

        self.gcn_hidden = gcn_hidden
        self.gcn1 = DenseGCNLayer(self.node_in, gcn_hidden)
        self.gcn2 = DenseGCNLayer(gcn_hidden, gcn_hidden)

        self.edge_mlp = nn.Sequential(
            nn.Linear(2 * gcn_hidden + edge_feat_dim, edge_hidden),
            nn.ReLU(),
            nn.Linear(edge_hidden, 1),
        )

    def encode_nodes(self, seq_tokens, x_cov):
        E = self.embed(seq_tokens)
        mask = (seq_tokens != 0).float().unsqueeze(-1)
        if self.use_lstm:
            lengths = mask.squeeze(-1).sum(dim=1).clamp(min=1).long()
            H0, _ = self.lstm(E)
            Hseq = (H0 * mask).sum(dim=1) / lengths.clamp(min=1).unsqueeze(-1).float()
        else:
            lengths = mask.squeeze(-1).sum(dim=1).clamp(min=1)
            Hseq = (E * mask).sum(dim=1) / lengths.unsqueeze(-1)
        X = torch.cat([Hseq, x_cov], dim=1)
        return X

    def forward(self, data: Data):
        device = data.seq_tokens.device
        N = data.num_nodes
        E = data.edge_index.size(1)
        X0 = self.encode_nodes(data.seq_tokens, data.x_cov)

        out_node = X0.new_zeros((N, self.gcn_hidden))
        batch_vec = data.batch if hasattr(data, "batch") else torch.zeros(N, dtype=torch.long, device=device)
        num_graphs = int(batch_vec.max().item()) + 1 if N > 0 else 0

        for g in range(num_graphs):
            node_ids = (batch_vec == g).nonzero(as_tuple=False).view(-1)
            n_nodes = int(node_ids.numel())
            if n_nodes == 0:
                continue

            local_map = torch.full((N,), -1, dtype=torch.long, device=device)
            local_map[node_ids] = torch.arange(n_nodes, device=device, dtype=torch.long)

            ei = data.edge_index
            keep = (local_map[ei[0]] >= 0) & (local_map[ei[1]] >= 0)
            keep_idx = keep.nonzero(as_tuple=False).view(-1)

            if keep_idx.numel() == 0:
                H = F.relu(self.gcn1(X0[node_ids], torch.empty(2,0,device=device,dtype=torch.long), n_nodes=n_nodes))
                H = F.relu(self.gcn2(H,          torch.empty(2,0,device=device,dtype=torch.long), n_nodes=n_nodes))
                out_node[node_ids] = H
                continue

            src_local = local_map[ei[0, keep_idx]]
            dst_local = local_map[ei[1, keep_idx]]
            edge_local = torch.stack([src_local, dst_local], dim=0)

            H = F.relu(self.gcn1(X0[node_ids], edge_local, n_nodes))
            H = F.relu(self.gcn2(H,            edge_local, n_nodes))
            out_node[node_ids] = H

        if E == 0:
            return torch.empty((0,), device=device)

        u, v = data.edge_index
        U = out_node[u]
        V = out_node[v]
        EA = data.edge_attr if hasattr(data, "edge_attr") and data.edge_attr.numel() > 0 \
             else out_node.new_zeros((E, 5))
        feats = torch.cat([U, V, EA], dim=1)
        logits = self.edge_mlp(feats).squeeze(-1)  # [E]
        return logits

# --------------------------
# Inference helpers
# --------------------------

@torch.no_grad()
def infer_and_emit(
    model: nn.Module,
    loader: DataLoader,
    device: torch.device,
    out_path: Path,
    keep_all_sinks: bool = True,
) -> None:
    """
    For each bubble graph:
      - score edges
      - for every source node with >=2 outgoing edges, keep only argmax(logit)
      - optionally keep edges from sources with a single outgoing edge (usually helpful)
      - write one JSON record per bubble:
          { "bubble_id": ..., "keep_edges": [ {"u_seq":..., "v_seq":...}, ... ] }
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        for batch in loader:
            batch = batch.to(device)
            logits = model(batch)
            if logits.numel() == 0:
                # still emit empty decision objects for consistency
                _emit_empty_decisions(batch, fh)
                continue

            u = batch.edge_index[0]
            v = batch.edge_index[1]
            node_batch = batch.batch if hasattr(batch, "batch") else torch.zeros(batch.num_nodes, dtype=torch.long, device=device)
            edge_batch = node_batch[u] if u.numel() else torch.zeros((0,), dtype=torch.long, device=device)
            num_graphs = int(node_batch.max().item()) + 1 if batch.num_nodes > 0 else 0

            # split by graph
            for g in range(num_graphs):
                g_edge_idx = (edge_batch == g).nonzero(as_tuple=False).view(-1)
                if g_edge_idx.numel() == 0:
                    rec = {
                        "bubble_id": _safe_get_attr(batch, "bubble_id", g),
                        "keep_edges": []
                    }
                    fh.write(json.dumps(rec) + "\n")
                    continue

                # For decision points: group edges by same source node
                src_nodes = u[g_edge_idx]
                dst_nodes = v[g_edge_idx]
                edge_scores = logits[g_edge_idx]

                # Build groups { source_node_idx -> [(edge_idx_in_g_edge_idx, score), ...] }
                groups: Dict[int, List[Tuple[int, float]]] = {}
                for i in range(g_edge_idx.numel()):
                    s = int(src_nodes[i].item())
                    groups.setdefault(s, []).append((i, float(edge_scores[i].item())))

                keep_pairs: List[Tuple[int, int]] = []
                for s, items in groups.items():
                    if len(items) >= 2:
                        # choose best
                        best_local_idx = max(items, key=lambda t: t[1])[0]
                        keep_pairs.append((s, int(dst_nodes[best_local_idx].item())))
                    elif keep_all_sinks and len(items) == 1:
                        # keep lone edge if requested
                        keep_pairs.append((s, int(dst_nodes[items[0][0]].item())))

                # map node indices back to sequences
                node_seqs = _node_seqs_for_graph(batch, g)
                keep_edges_json = []
                for (s_idx, t_idx) in keep_pairs:
                    try:
                        u_seq = node_seqs[s_idx]
                        v_seq = node_seqs[t_idx]
                    except Exception:
                        # fallback: skip if indices out of range
                        continue
                    keep_edges_json.append({"u_seq": u_seq, "v_seq": v_seq})

                rec = {
                    "bubble_id": _safe_get_attr(batch, "bubble_id", g),
                    "keep_edges": keep_edges_json
                }
                fh.write(json.dumps(rec) + "\n")

def _safe_get_attr(batch: Data, name: str, default_val):
    try:
        val = getattr(batch, name)
        # If batched tensors, we can't reliably split without extra bookkeeping; use default.
        if isinstance(val, torch.Tensor):
            return default_val
        return val
    except Exception:
        return default_val

def _node_seqs_for_graph(batch: Data, g: int) -> List[str]:

    device = batch.seq_tokens.device
    node_batch = batch.batch if hasattr(batch, "batch") else torch.zeros(batch.num_nodes, dtype=torch.long, device=device)

    node_ids = (node_batch == g).nonzero(as_tuple=False).view(-1)

    ns = batch.node_seqs

    if not ns:
        return []

    if isinstance(ns[0], str):
        flat = ns
    elif isinstance(ns[0], list):
        flat = []
        for lst in ns:
            flat.extend(lst)
    else:
        return []

    out: List[str] = []
    for i in range(int(node_ids.numel())):
        gi = int(node_ids[i].item())
        if gi < 0 or gi >= len(flat):
            continue
        out.append(flat[gi])
    return out


def _emit_empty_decisions(batch: Data, fh):
    device = batch.seq_tokens.device
    node_batch = batch.batch if hasattr(batch, "batch") else torch.zeros(batch.num_nodes, dtype=torch.long, device=device)
    num_graphs = int(node_batch.max().item()) + 1 if batch.num_nodes > 0 else 0
    for g in range(num_graphs):
        rec = {"bubble_id": _safe_get_attr(batch, "bubble_id", g), "keep_edges": []}
        fh.write(json.dumps(rec) + "\n")

# --------------------------
# CLI / main
# --------------------------

def load_device(use_directml: bool, force_cpu: bool) -> torch.device:
    if force_cpu:
        return torch.device("cpu")
    if use_directml:
        try:
            import torch_directml
            return torch_directml.device()
        except Exception:
            pass
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")

def load_model(ckpt_path: Path, device: torch.device, use_lstm: bool = False) -> nn.Module:
    model = HyperbubbleGNN(
        vocab_size=5,
        embed_dim=16,
        gcn_hidden=64,
        edge_hidden=64,
        edge_feat_dim=5,
        use_lstm=use_lstm,
    ).to(device)

    # Accept both formats:
    #   1) {"state_dict": <...>, ...}
    #   2) raw state_dict
    ckpt = torch.load(ckpt_path, map_location="cpu")
    if isinstance(ckpt, dict) and "state_dict" in ckpt:
        state_dict = ckpt["state_dict"]
    else:
        state_dict = ckpt
    model.load_state_dict(state_dict, strict=False)
    model.eval()
    return model

def main():
    ap = argparse.ArgumentParser(description="Run Hyperbubble GNN on BubbleGun JSONL and emit decisions JSONL.")
    ap.add_argument("--in", dest="inp", required=True, help="Input bubbles JSONL (from BubbleGun).")
    ap.add_argument("--out", dest="out", required=True, help="Output decisions JSONL for Rust resolver.")
    ap.add_argument("--ckpt", required=True, help="Path to saved model checkpoint/state_dict.")
    ap.add_argument("--batch-size", type=int, default=1, help="Inference batch size.")
    ap.add_argument("--num-workers", type=int, default=0, help="DataLoader workers.")
    ap.add_argument("--directml", action="store_true", help="Use DirectML if available.")
    ap.add_argument("--cpu", action="store_true", help="Force CPU.")
    ap.add_argument("--use-lstm", action="store_true", help="Use LSTM token encoder (must match training).")
    ap.add_argument("--gfa", required=True, help="GFA graph used with BubbleGun JSON to map node IDs to k-mers.")
    args = ap.parse_args()

    device = load_device(use_directml=args.directml, force_cpu=args.cpu)
    model = load_model(Path(args.ckpt), device=device, use_lstm=args.use_lstm)

    # Build graphs directly from BubbleGun JSON + GFA
    ds = BubbleGunBubbleDataset(args.inp, args.gfa)
    loader = DataLoader(ds, batch_size=args.batch_size, shuffle=False, num_workers=args.num_workers)

    out_path = Path(args.out)
    infer_and_emit(model, loader, device, out_path)
    print(f"[gnn-infer] wrote decisions to {out_path}")

if __name__ == "__main__":
    main()
