#!/usr/bin/env python3
import argparse
import json
import math
from typing import Any, Dict, Iterable, List, Optional, Tuple

# ---------- IO ----------
def iter_jsonl(path: str) -> Iterable[Dict[str, Any]]:
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            yield json.loads(line)

# ---------- Coverage helpers ----------
def edge_cov(e: Dict[str, Any]) -> Optional[float]:
    """Extract a single coverage value per edge (prefers cov_mean; falls back to cov_min)."""
    if "cov_mean" in e and e["cov_mean"] is not None:
        return float(e["cov_mean"])
    if "cov_min" in e and e["cov_min"] is not None:
        return float(e["cov_min"])
    return None

# ---------- Decision mining (orientation-agnostic) ----------
def decisions_for_record(rec: Dict[str, Any]) -> List[Tuple[int, List[int], int, Optional[float]]]:
    """
    Return a list of decisions for this bubble.
    Each item is: (source_node_key, cand_edge_indices, gold_pos, ratio)
      - source_node_key is an integer index into a deduplicated source list (not used for scoring, informational)
      - cand_edge_indices: indices of all outgoing edges from that decision node
      - gold_pos: index *within* cand_edge_indices that is labeled (supervision)
      - ratio: local difficulty ratio = max(labeled_cov, best_other_cov) / min(labeled_cov, best_other_cov)
               or None if coverage missing or non-positive
    A node is a decision node iff it has >=2 outgoing edges and exactly one of them is labeled.
    """
    edges: List[Dict[str, Any]] = rec.get("edges", [])
    paths: Any = rec.get("paths", [])
    lp = rec.get("label_path")
    if not isinstance(paths, list) or lp is None or not (0 <= int(lp) < len(paths)):
        return []

    labeled_edge_set = set(paths[int(lp)])

    # Build outgoing adjacency by source_seq (strings)
    from collections import defaultdict
    out_by_src = defaultdict(list)  # src_seq -> list of edge indices
    for i, e in enumerate(edges):
        out_by_src[e["source_seq"]].append(i)

    # Source nodes that are on the labeled path as *sources*
    labeled_sources = set(edges[i]["source_seq"] for i in labeled_edge_set if 0 <= i < len(edges))

    # Freeze order for reproducibility
    labeled_sources = list(labeled_sources)

    decisions: List[Tuple[int, List[int], int, Optional[float]]] = []
    for si, src_seq in enumerate(labeled_sources):
        cand = out_by_src.get(src_seq, [])
        if len(cand) < 2:
            continue  # not a branching node

        # Which candidates are labeled?
        labeled_mask = [int(ci in labeled_edge_set) for ci in cand]
        if sum(labeled_mask) != 1:
            # ambiguous supervision -> skip
            continue

        gold_pos = labeled_mask.index(1)

        # Compute local ratio for difficulty analysis
        # Use only edges with valid positive coverage
        covs = []
        for ci in cand:
            c = edge_cov(edges[ci])
            covs.append(c if (c is not None) else None)

        labeled_cov = covs[gold_pos]
        # best "other" is the max coverage among non-gold candidates
        other_covs = [covs[j] for j in range(len(cand)) if j != gold_pos]
        best_other_cov = None
        if other_covs:
            # best among valid coverage values
            filtered = [x for x in other_covs if x is not None]
            best_other_cov = max(filtered) if filtered else None

        ratio = None
        if labeled_cov is not None and best_other_cov is not None and labeled_cov > 0 and best_other_cov > 0:
            big = max(labeled_cov, best_other_cov)
            small = min(labeled_cov, best_other_cov)
            ratio = big / small if small > 0 else None

        decisions.append((si, cand, gold_pos, ratio))

    return decisions

# ---------- Simple coverage chooser (local) ----------
def predict_local_by_cov(edges: List[Dict[str, Any]], cand: List[int]) -> Optional[int]:
    """
    Predict the best candidate index by maximum coverage among 'cand'.
    Returns the index within 'cand' or None if coverage is missing.
    """
    best_idx = None
    best_val = None
    for j, ei in enumerate(cand):
        c = edge_cov(edges[ei])
        if c is None:
            continue
        if (best_val is None) or (c > best_val):
            best_val = c
            best_idx = j
    return best_idx

# ---------- Loading and per-bubble expansion ----------
def load_expand_decisions(paths: List[str]) -> List[Dict[str, Any]]:
    """
    Expand bubbles into decision records.
    Each returned dict has:
        - 'bubble_id', 'k', ... (passthrough if present)
        - '_cand': List[int] of candidate edge indices at this decision
        - '_gold_pos': int in [0, len(_cand)-1]
        - '_pred_pos': Optional[int]
        - '_correct': Optional[bool]
        - '_ratio': Optional[float]
        - '_bubble_ok': Optional[bool] (filled later)
        - '_source_seq': str (source node sequence for this decision)
    """
    out: List[Dict[str, Any]] = []
    for p in paths:
        for rec in iter_jsonl(p):
            decs = decisions_for_record(rec)
            if not decs:
                continue
            edges = rec.get("edges", [])
            bubble_ok_all = True  # temporary; we compute after pred
            for si, cand, gold_pos, ratio in decs:
                pred_pos = predict_local_by_cov(edges, cand)
                correct = (pred_pos == gold_pos) if (pred_pos is not None) else None

                # find the source node sequence for logging
                src_seq = edges[cand[0]]["source_seq"] if cand else None

                item = dict(rec)  # shallow copy; OK for read-only
                item["_cand"] = cand
                item["_gold_pos"] = gold_pos
                item["_pred_pos"] = pred_pos
                item["_correct"] = correct
                item["_ratio"] = ratio
                item["_source_seq"] = src_seq
                out.append(item)
    # Compute bubble-level correctness: a bubble is correct iff all its decisions are correct (not None)
    # Group by bubble_id if present else by object id
    from collections import defaultdict
    group = defaultdict(list)
    for i, r in enumerate(out):
        key = r.get("bubble_id", f"idx:{i}")
        group[key].append(i)
    for key, idxs in group.items():
        # If any decision has _correct == False, bubble is incorrect.
        # If any decision has _correct is None, we mark bubble as None (unknown).
        bubble_ok = True
        unknown = False
        for i in idxs:
            c = out[i]["_correct"]
            if c is None:
                unknown = True
            elif c is False:
                bubble_ok = False
        for i in idxs:
            out[i]["_bubble_ok"] = (None if unknown else bubble_ok)
    return out

# ---------- Threshold scan ----------
def frange(start: float, stop: float, step: float) -> List[float]:
    vals: List[float] = []
    x = start
    while x <= stop + 1e-12:
        vals.append(float(round(x, 6)))
        x += step
    return vals

def evaluate_below_threshold_decisions(
    records: List[Dict[str, Any]],
    thresholds: List[float],
    nontrivial_cap: float
) -> List[Dict[str, Any]]:
    """
    Decision-level accuracy below threshold.
    We keep decisions where _ratio is known and <= nontrivial_cap (difficulty ceiling),
    then for each threshold tau compute accuracy over subset with _ratio <= tau.
    """
    results: List[Dict[str, Any]] = []
    pool = [r for r in records if r.get("_ratio") is not None and r["_ratio"] <= nontrivial_cap]
    missing_ratio_count = sum(1 for r in records if r.get("_ratio") is None)

    for tau in thresholds:
        subset = [r for r in pool if r["_ratio"] <= tau]
        denom = sum(1 for r in subset if r.get("_correct") is not None)
        num_ok = sum(1 for r in subset if r.get("_correct") is True)
        acc = (num_ok / denom) if denom > 0 else float("nan")
        results.append({
            "threshold": tau,
            "count": len(subset),          # number of decisions considered
            "eval_count": denom,           # decisions with known prediction
            "accuracy": acc,
            "missing_ratio_global": missing_ratio_count
        })
    return results

def bubble_level_summary(records: List[Dict[str, Any]]) -> Tuple[int, int, int, float]:
    """
    Summarize bubble-level metrics from decision-expanded records.
    Returns: (num_bubbles_with_decisions, num_bubbles_eval, num_bubbles_correct, bubble_acc)
    """
    from collections import defaultdict
    bubbles = defaultdict(list)
    for i, r in enumerate(records):
        key = r.get("bubble_id", f"idx:{i}")
        bubbles[key].append(i)

    with_decisions = len(bubbles)
    evaluable = 0
    correct = 0
    for key, idxs in bubbles.items():
        # a bubble is evaluable if none of its decisions have _correct is None
        any_unknown = any(records[i]["_correct"] is None for i in idxs)
        if any_unknown:
            continue
        evaluable += 1
        all_ok = all(records[i]["_correct"] is True for i in idxs)
        if all_ok:
            correct += 1
    acc = (correct / evaluable) if evaluable > 0 else float("nan")
    return with_decisions, evaluable, correct, acc

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(description="Decision-level coverage heuristic evaluation (all decision nodes).")
    ap.add_argument("files", nargs="+", help="Input JSONL files (bubbles).")
    ap.add_argument("--start", type=float, default=1.0, help="Start threshold (default: 1.0).")
    ap.add_argument("--stop", type=float, default=21.0, help="Stop threshold (default: 21.0).")
    ap.add_argument("--step", type=float, default=0.5, help="Step size for thresholds.")
    ap.add_argument("--nontrivial-cap", type=float, default=21.0, help="Upper cap for difficulty (default: 21.0).")
    args = ap.parse_args()

    # Expand bubbles into per-decision records
    decision_records = load_expand_decisions(args.files)

    # Decision-level threshold curve
    thresholds = frange(args.start, args.stop, args.step)
    results = evaluate_below_threshold_decisions(decision_records, thresholds, args.nontrivial_cap)

    # Bubble-level summary
    bubbles_with_decisions, bubbles_eval, bubbles_ok, bub_acc = bubble_level_summary(decision_records)

    print(f"# decision-level coverage accuracy (ratio ≤ t), non-trivial cap ≤ {args.nontrivial_cap:g}")
    header = f"{'threshold':>9}  {'count':>7}  {'eval_count':>11}  {'accuracy':>9}"
    print(header)
    print("-" * len(header))
    for r in results:
        acc = "nan" if math.isnan(r["accuracy"]) else f"{100.0 * r['accuracy']:.2f}%"
        print(f"{r['threshold']:9.1f}  {r['count']:7d}  {r['eval_count']:11d}  {acc:>9}")

    print("\n# bubble-level (all decision nodes must be correct)")
    if math.isnan(bub_acc):
        print(f"bubbles_with_decisions={bubbles_with_decisions}  evaluable={bubbles_eval}  correct={bubbles_ok}  acc=nan")
    else:
        print(f"bubbles_with_decisions={bubbles_with_decisions}  evaluable={bubbles_eval}  correct={bubbles_ok}  acc={bub_acc:.3f}")

if __name__ == "__main__":
    main()
