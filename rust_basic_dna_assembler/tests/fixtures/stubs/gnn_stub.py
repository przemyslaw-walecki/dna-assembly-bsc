#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple, Any, Iterable, Optional


def parse_gfa(gfa_path: str) -> Tuple[Dict[str, str], Dict[Tuple[str, str], int], Dict[str, List[Tuple[str, int]]]]:
    """
    Reads a minimal subset of my GFA:
      S <nid> <seq> ... KC:i:<int>
      L <u> <+> <v> <+> <overlap> ... EC:i:<int>

    Returns:
      node_seq[nid] = seq
      edge_ec[(u,v)] = ec
      adj[u] = [(v, ec), ...]
    """
    node_seq: Dict[str, str] = {}
    edge_ec: Dict[Tuple[str, str], int] = {}
    adj: Dict[str, List[Tuple[str, int]]] = {}

    with open(gfa_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("H"):
                continue
            cols = line.split("\t")
            tag = cols[0]
            if tag == "S" and len(cols) >= 3:
                nid = cols[1]
                seq = cols[2]
                node_seq[nid] = seq
            elif tag == "L" and len(cols) >= 7:
                u = cols[1]
                v = cols[3]
                ec = 0
                for field in cols[6:]:
                    if field.lower().startswith("ec:i:"):
                        try:
                            ec = int(field.split(":")[-1])
                        except Exception:
                            ec = 0
                edge_ec[(u, v)] = ec
                adj.setdefault(u, []).append((v, ec))

    return node_seq, edge_ec, adj


def iter_bubbles(bubbles_path: str) -> Iterable[Dict[str, Any]]:
    """
    Matches the tolerant parsing pattern used in your current gnn_infer:
    - line can be: { "bubbles":[...] } OR dict containing dict-values OR list of dicts
    - each bubble is dict with:
        ends: [start_id, end_id]
        inside: [ ... ]
        id: optional
    Yields normalized records:
      { bubble_id, start, end, inside[list[str]] }
    """
    with open(bubbles_path, "r", encoding="utf-8") as fh:
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
                if not isinstance(bs, list):
                    continue
                for b in bs:
                    if not isinstance(b, dict):
                        continue
                    ends = b.get("ends", [])
                    inside = b.get("inside", [])
                    if not (isinstance(ends, list) and len(ends) == 2):
                        continue
                    if not isinstance(inside, list):
                        inside = []
                    bubble_id = b.get("id", None)
                    if bubble_id is None:
                        bubble_id = b.get("bubble_id", None)
                    if bubble_id is None:
                        bubble_id = 0  # will be overwritten by running index if needed
                    yield {
                        "bubble_id": bubble_id,
                        "start": str(ends[0]),
                        "end": str(ends[1]),
                        "inside": [str(x) for x in inside],
                    }


def choose_keep_edges_for_bubble(
    node_seq: Dict[str, str],
    adj: Dict[str, List[Tuple[str, int]]],
    bubble: Dict[str, Any],
) -> List[Tuple[str, str]]:
    """
    Deterministic decision rule:
      - Node set = inside + {start,end}
      - Consider all edges u->v with u,v in node set
      - Group by u
        - if >=2 outgoing in group: keep max EC; tie: min v_seq lexicographically
        - if ==1: keep that one
    Returns list of pairs (u_seq, v_seq).
    """
    node_ids = set(bubble["inside"])
    node_ids.add(bubble["start"])
    node_ids.add(bubble["end"])

    # Restrict adjacency to subgraph
    groups: Dict[str, List[Tuple[str, int]]] = {}
    for u in node_ids:
        for (v, ec) in adj.get(u, []):
            if v in node_ids:
                groups.setdefault(u, []).append((v, ec))

    keep: List[Tuple[str, str]] = []
    for u, outs in groups.items():
        if len(outs) == 0:
            continue
        if len(outs) == 1:
            v, _ = outs[0]
            u_seq = node_seq.get(u)
            v_seq = node_seq.get(v)
            if u_seq is not None and v_seq is not None:
                keep.append((u_seq, v_seq))
            continue

        # pick best by ec desc, then v_seq asc (deterministic tie-break)
        def key(item: Tuple[str, int]) -> Tuple[int, str]:
            v, ec = item
            v_seq = node_seq.get(v, "")
            return (ec, "\x00" if v_seq == "" else v_seq)

        best_ec = max(ec for (_, ec) in outs)
        best_candidates = [x for x in outs if x[1] == best_ec]
        best_candidates.sort(key=lambda x: node_seq.get(x[0], ""))  # v_seq asc
        v_best, _ = best_candidates[0]

        u_seq = node_seq.get(u)
        v_seq = node_seq.get(v_best)
        if u_seq is not None and v_seq is not None:
            keep.append((u_seq, v_seq))

    return keep


def main() -> None:
    ap = argparse.ArgumentParser(description="Deterministic stub for gnn_infer.py (no ML).")
    ap.add_argument("--in", dest="inp", required=True, help="Input bubbles JSONL (from BubbleGun).")
    ap.add_argument("--out", dest="out", required=True, help="Output decisions JSONL for Rust resolver.")
    ap.add_argument("--gfa", required=True, help="GFA graph used with BubbleGun JSON to map node IDs to k-mers.")
    # Ignored but accepted for CLI compatibility:
    ap.add_argument("--ckpt", required=False, default="", help="Ignored (stub).")
    ap.add_argument("--batch-size", type=int, default=1, help="Ignored (stub).")
    ap.add_argument("--num-workers", type=int, default=0, help="Ignored (stub).")
    ap.add_argument("--directml", action="store_true", help="Ignored (stub).")
    ap.add_argument("--cpu", action="store_true", help="Ignored (stub).")
    args = ap.parse_args()

    node_seq, _edge_ec, adj = parse_gfa(args.gfa)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Ensure deterministic bubble_id assignment if input lacks stable ids:
    # keep running counter and overwrite bubble_id when it was default-ish.
    bubble_idx = 0

    with out_path.open("w", encoding="utf-8") as fh:
        for bub in iter_bubbles(args.inp):
            # normalize bubble_id
            bid = bub.get("bubble_id", 0)
            if not isinstance(bid, int):
                try:
                    bid = int(bid)
                except Exception:
                    bid = bubble_idx
            if bid == 0 and bubble_idx != 0:
                bid = bubble_idx

            keep_pairs = choose_keep_edges_for_bubble(node_seq, adj, bub)
            keep_edges_json = [{"u_seq": u, "v_seq": v} for (u, v) in keep_pairs]

            rec = {"bubble_id": bid, "keep_edges": keep_edges_json}
            fh.write(json.dumps(rec) + "\n")

            bubble_idx += 1

    print(f"[gnn-stub] wrote decisions to {out_path}")


if __name__ == "__main__":
    main()
