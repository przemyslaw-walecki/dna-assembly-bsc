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
from torch_geometric.nn import GCNConv

# --------------------------
# Utilities / tokenization
# --------------------------

NUC2ID = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 5}
ID2NUC = {v: k for k, v in NUC2ID.items()}

def seq_to_tokens(seq: str) -> torch.LongTensor:
    # NOTE: 0 is reserved as "padding", but here we don't pad; values are in {1..5}
    return torch.tensor([NUC2ID.get(c, 5) for c in seq], dtype=torch.long)

def _safe_get(d: Dict[str, Any], key: str, default):
    v = d.get(key, default)
    return default if v is None else v

def tokens_to_onehot(tokens: torch.LongTensor, num_classes: int = 6) -> torch.Tensor:
    """
    DML-friendly one-hot via broadcasting.
    tokens: [N, K] long
    returns: [N, K, num_classes] float
    """
    if tokens.numel() == 0:
        # preserve shape for empty graphs
        return tokens.new_zeros((*tokens.shape, num_classes), dtype=torch.float32)

    classes = torch.arange(num_classes, device=tokens.device, dtype=tokens.dtype).view(1, 1, num_classes)
    return (tokens.unsqueeze(-1) == classes).to(torch.float32)

# --------------------------
# Dataset for BubbleGun JSONL
# --------------------------

class BubbleGunBubbleDataset(Dataset):
    """
    Dataset built directly from:
      - BubbleGun JSON/JSONL (with numeric node IDs)
      - Graph GFA (from Rust assembler) for seq + coverage + adjacency

    For each bubble:
      - nodes = inside + endpoints
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
                continue
            local_idx = len(node_seqs)
            id_to_local[nid] = local_idx
            node_seqs.append(seq)
            cov = self.node_cov.get(nid, 0.0)
            node_covs.append(float(cov))

        n_nodes = len(node_seqs)
        if n_nodes == 0:
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
                edge_src.append(u_idx)
                edge_dst.append(v_idx)
                edge_attr.append([1.0, float(self.k), ec, ec, 0.0])

        if edge_src:
            edge_index = torch.tensor([edge_src, edge_dst], dtype=torch.long)
            edge_attr_t = torch.tensor(edge_attr, dtype=torch.float32)
        else:
            edge_index = torch.zeros((2, 0), dtype=torch.long)
            edge_attr_t = torch.zeros((0, 5), dtype=torch.float32)

        seq_tokens = torch.stack([seq_to_tokens(s) for s in node_seqs], dim=0)  # [N,K]
        x_cov = torch.tensor(node_covs, dtype=torch.float32).unsqueeze(1)       # [N,1]

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
# Model (NEW): CNN(one-hot) + 2xGCNConv + edge MLP
# --------------------------

class HyperbubbleGNN(nn.Module):
    """
    CNN-based sequence encoder + 2× GCNConv + MLP edge classifier.
    Matches your new architecture (one-hot -> Conv1d -> pooling -> concat cov).
    """
    def __init__(
        self,
        vocab_size: int = 6,     # 0..5, where 1..5 used for A,C,G,T,N
        cnn_channels: int = 32,
        gcn_hidden: int = 64,
        edge_hidden: int = 64,
        edge_feat_dim: int = 5,
        dropout: float = 0.0,
    ):
        super().__init__()

        self.vocab_size = vocab_size
        self.cnn_channels = cnn_channels
        self.edge_feat_dim = edge_feat_dim

        self.cnn = nn.Conv1d(vocab_size, cnn_channels, kernel_size=3, padding=1)
        self.node_in = cnn_channels + 1

        self.gcn1 = GCNConv(self.node_in, gcn_hidden)
        self.gcn2 = GCNConv(gcn_hidden, gcn_hidden)

        self.dropout = nn.Dropout(dropout)

        self.edge_mlp = nn.Sequential(
            nn.Linear(2 * gcn_hidden + edge_feat_dim, edge_hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(edge_hidden, 1),
        )

    def encode_nodes(self, seq_tokens: torch.Tensor, x_cov: torch.Tensor) -> torch.Tensor:
        """
        seq_tokens: [N,K] long
        x_cov:      [N,1] float
        returns:    [N, cnn_channels+1]
        """
        if seq_tokens.numel() == 0:
            return x_cov.new_zeros((0, self.cnn_channels + 1))

        onehot = tokens_to_onehot(seq_tokens, num_classes=self.vocab_size)  # [N,K,V]
        x = onehot.permute(0, 2, 1)                                        # [N,V,K]

        h = F.relu(self.cnn(x))                                            # [N,C,K]
        h = h.mean(dim=2)                                                  # [N,C]
        return torch.cat([h, x_cov], dim=1)                                # [N,C+1]

    def forward(self, data: Data) -> torch.Tensor:
        """
        returns logits [E] (one per edge in data.edge_index)
        """
        X0 = self.encode_nodes(data.seq_tokens, data.x_cov)  # [N,node_in]

        if data.num_nodes == 0:
            return X0.new_zeros((0,))

        # Node propagation on sparse edge_index
        H = F.relu(self.gcn1(X0, data.edge_index))
        H = self.dropout(H)
        H = F.relu(self.gcn2(H, data.edge_index))
        H = self.dropout(H)

        E = data.edge_index.size(1)
        if E == 0:
            return H.new_zeros((0,))

        u, v = data.edge_index
        U = H[u]
        V = H[v]

        if hasattr(data, "edge_attr") and data.edge_attr.numel() > 0:
            EA = data.edge_attr
        else:
            EA = H.new_zeros((E, self.edge_feat_dim))

        feats = torch.cat([U, V, EA], dim=1)                 # [E, 2*gcn_hidden + edge_feat_dim]
        logits = self.edge_mlp(feats).squeeze(-1)            # [E]
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
      - optionally keep edges from sources with a single outgoing edge
      - write one JSON record per bubble:
          { "bubble_id": ..., "keep_edges": [ {"u_seq":..., "v_seq":...}, ... ] }
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        for batch in loader:
            batch = batch.to(device)
            logits = model(batch)

            if logits.numel() == 0:
                _emit_empty_decisions(batch, fh)
                continue

            u = batch.edge_index[0]
            v = batch.edge_index[1]
            node_batch = batch.batch if hasattr(batch, "batch") else torch.zeros(batch.num_nodes, dtype=torch.long, device=device)
            edge_batch = node_batch[u] if u.numel() else torch.zeros((0,), dtype=torch.long, device=device)
            num_graphs = int(node_batch.max().item()) + 1 if batch.num_nodes > 0 else 0

            for g in range(num_graphs):
                g_edge_idx = (edge_batch == g).nonzero(as_tuple=False).view(-1)
                if g_edge_idx.numel() == 0:
                    rec = {"bubble_id": _safe_get_attr(batch, "bubble_id", g), "keep_edges": []}
                    fh.write(json.dumps(rec) + "\n")
                    continue

                src_nodes = u[g_edge_idx]
                dst_nodes = v[g_edge_idx]
                edge_scores = logits[g_edge_idx]

                groups: Dict[int, List[Tuple[int, float]]] = {}
                for i in range(g_edge_idx.numel()):
                    s = int(src_nodes[i].item())
                    groups.setdefault(s, []).append((i, float(edge_scores[i].item())))

                keep_pairs: List[Tuple[int, int]] = []
                for s, items in groups.items():
                    if len(items) >= 2:
                        best_local_idx = max(items, key=lambda t: t[1])[0]
                        keep_pairs.append((s, int(dst_nodes[best_local_idx].item())))
                    elif keep_all_sinks and len(items) == 1:
                        keep_pairs.append((s, int(dst_nodes[items[0][0]].item())))

                node_seqs = _node_seqs_for_graph(batch, g)
                keep_edges_json = []
                for (s_idx, t_idx) in keep_pairs:
                    try:
                        u_seq = node_seqs[s_idx]
                        v_seq = node_seqs[t_idx]
                    except Exception:
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

def load_model(ckpt_path: Path, device: torch.device) -> nn.Module:
    model = HyperbubbleGNN(
        vocab_size=6,
        cnn_channels=8,
        gcn_hidden=16,
        edge_hidden=256,
        edge_feat_dim=5,
        dropout=0.0,
    ).to(device)

    ckpt = torch.load(ckpt_path, map_location="cpu")
    state_dict = ckpt["state_dict"] if isinstance(ckpt, dict) and "state_dict" in ckpt else ckpt
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
    ap.add_argument("--gfa", required=True, help="GFA graph used with BubbleGun JSON to map node IDs to k-mers.")
    args = ap.parse_args()

    device = load_device(use_directml=args.directml, force_cpu=args.cpu)
    model = load_model(Path(args.ckpt), device=device)

    ds = BubbleGunBubbleDataset(args.inp, args.gfa)
    loader = DataLoader(ds, batch_size=args.batch_size, shuffle=False, num_workers=args.num_workers)

    out_path = Path(args.out)
    infer_and_emit(model, loader, device, out_path)
    print(f"[gnn-infer] wrote decisions to {out_path}")

if __name__ == "__main__":
    main()
