#!/usr/bin/env python3
"""
Improved heuristic evaluator for bubble decisions.

Changes vs original:
 - cleaner structure
 - accepts multiple input files
 - supports per-genome and global aggregation
 - identical scoring logic
"""

import argparse
import json
import math
from typing import Dict, Any, List, Iterable, Tuple, Optional
from collections import defaultdict

# ============================================================
# IO
# ============================================================

def iter_jsonl(path: str) -> Iterable[Dict[str, Any]]:
    """Iterate over JSONL records."""
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                yield json.loads(line)

# ============================================================
# Coverage helpers
# ============================================================

def edge_cov(e: Dict[str, Any]) -> Optional[float]:
    """Extract coverage for an edge, preferring cov_mean."""
    if e.get("cov_mean") is not None:
        return float(e["cov_mean"])
    if e.get("cov_min") is not None:
        return float(e["cov_min"])
    return None


def average_cov_for_path(edges: List[Dict[str, Any]], path: List[int]) -> Optional[float]:
    """Average coverage for a path: prefer cov_mean, fallback to cov_min for each edge."""
    covs = []
    for idx in path:
        if not (0 <= idx < len(edges)):
            return None
        e = edges[idx]
        cov = e.get("cov_mean", e.get("cov_min", None))
        if cov is None:
            return None
        covs.append(float(cov))
    if not covs:
        return None
    return sum(covs) / len(covs)


def coverage_ratio_two_paths(edges: List[Dict[str, Any]], paths: List[List[int]]) -> Optional[float]:
    """
    Compute ratio (≥1.0) = max(path_avg) / min(path_avg) for exactly two paths.
    Returns None if not computable or invalid.

    This is designed to match the logic used in the extractor script.
    """
    if not isinstance(paths, list) or len(paths) != 2:
        return None
    a = average_cov_for_path(edges, paths[0])
    b = average_cov_for_path(edges, paths[1])
    if a is None or b is None or a <= 0 or b <= 0:
        return None
    big, small = (a, b) if a >= b else (b, a)
    return big / small

# ============================================================
# Decision discovery
# ============================================================

def decisions_for_record(rec: Dict[str, Any]) -> List[Tuple[int, List[int], int, Optional[float]]]:
    """
    Return decisions for a single bubble record.

    Output: list of tuples:
      (source_key, cand_edge_indices, gold_pos, ratio)

    ratio is now the PATH-LEVEL coverage ratio between the two paths in this bubble,
    identical for all decisions originating from this bubble.
    """
    edges = rec.get("edges", [])
    paths = rec.get("paths", [])
    lp = rec.get("label_path")
    if not isinstance(paths, list) or lp is None:
        return []

    lp = int(lp)
    if lp < 0 or lp >= len(paths):
        return []

    labeled_edges = set(paths[lp])

    # compute a single path-level ratio for the whole bubble
    path_ratio = coverage_ratio_two_paths(edges, paths)

    out_by_src = defaultdict(list)
    for idx, e in enumerate(edges):
        out_by_src[e["source_seq"]].append(idx)

    labeled_sources = set(edges[i]["source_seq"] for i in labeled_edges if 0 <= i < len(edges))
    labeled_sources = list(labeled_sources)

    decisions = []
    for si, src in enumerate(labeled_sources):
        cand = out_by_src.get(src, [])
        if len(cand) < 2:
            continue

        labeled_mask = [1 if c in labeled_edges else 0 for c in cand]
        if sum(labeled_mask) != 1:
            continue

        gold_pos = labeled_mask.index(1)

        # ratio is the same for all decisions in this bubble:
        ratio = path_ratio

        decisions.append((si, cand, gold_pos, ratio))

    return decisions

# ============================================================
# Heuristic predictor
# ============================================================

def predict_local_by_cov(edges: List[Dict[str, Any]], cand: List[int]) -> Optional[int]:
    """
    Simple heuristic: pick the candidate edge with the highest coverage.
    Coverage is still per-edge here; only the "difficulty ratio" is path-level.
    """
    best_j = None
    best_v = None
    for j, ei in enumerate(cand):
        v = edge_cov(edges[ei])
        if v is None:
            continue
        if best_v is None or v > best_v:
            best_v = v
            best_j = j
    return best_j

# ============================================================
# Bubble expansion to decision-level records
# ============================================================

def load_expand_decisions(paths: List[str]) -> List[Dict[str, Any]]:
    """
    Expand all bubbles → decision records.
    """
    out = []

    for p in paths:
        for rec in iter_jsonl(p):
            decs = decisions_for_record(rec)
            edges = rec.get("edges", [])
            if not decs:
                continue
            for si, cand, gold_pos, ratio in decs:
                pred_pos = predict_local_by_cov(edges, cand)
                correct = None if pred_pos is None else (pred_pos == gold_pos)

                src_seq = edges[cand[0]]["source_seq"] if cand else None
                item = {
                    "bubble_id": rec.get("bubble_id"),
                    "_cand": cand,
                    "_gold_pos": gold_pos,
                    "_pred_pos": pred_pos,
                    "_correct": correct,
                    "_ratio": ratio,
                    "_source_seq": src_seq,
                }
                out.append(item)

    # compute bubble-level correctness
    grouped = defaultdict(list)
    for i, r in enumerate(out):
        key = r.get("bubble_id", f"idx:{i}")
        grouped[key].append(i)

    for key, idxs in grouped.items():
        unknown = any(out[i]["_correct"] is None for i in idxs)
        if unknown:
            bubble_val = None
        else:
            bubble_val = all(out[i]["_correct"] is True for i in idxs)
        for i in idxs:
            out[i]["_bubble_ok"] = bubble_val

    return out

# ============================================================
# Metrics
# ============================================================

def frange(start, stop, step):
    vals = []
    x = start
    while x <= stop + 1e-12:
        vals.append(round(x, 6))
        x += step
    return vals

def evaluate_thresholds(records, thresholds, cap):
    # Only consider records with a defined ratio and within the cap
    pool = [r for r in records if r["_ratio"] is not None and r["_ratio"] <= cap]
    global_missing = sum(1 for r in records if r["_ratio"] is None)

    results = []
    for t in thresholds:
        subset = [r for r in pool if r["_ratio"] <= t]
        eval_count = sum(1 for r in subset if r["_correct"] is not None)
        num_ok = sum(1 for r in subset if r["_correct"] is True)
        acc = num_ok / eval_count if eval_count > 0 else float("nan")
        results.append({
            "threshold": t,
            "count": len(subset),
            "eval_count": eval_count,
            "accuracy": acc,
            "missing_ratio": global_missing,
        })
    return results

def bubble_level(records):
    grouped = defaultdict(list)
    for i, r in enumerate(records):
        grouped[r["bubble_id"]].append(i)

    with_decs = len(grouped)
    evaluable = 0
    correct = 0
    for key, idxs in grouped.items():
        if any(records[i]["_correct"] is None for i in idxs):
            continue
        evaluable += 1
        if all(records[i]["_correct"] is True for i in idxs):
            correct += 1

    acc = correct / evaluable if evaluable > 0 else float("nan")
    return with_decs, evaluable, correct, acc

# ============================================================
# Genome aggregation
# ============================================================

def genome_id_from_path(path: str) -> str:
    """Extract genome name from filename like: ecoli_dataset_cov20_ratio_lt.jsonl"""
    name = path.split("/")[-1]
    return name.split("_dataset_cov")[0]

def group_records_by_genome(records: List[Dict[str, Any]], paths: List[str]):
    genome_to_records = defaultdict(list)
    for rec, p in zip(records, paths):
        gid = genome_id_from_path(p)
        genome_to_records[gid].append(rec)
    return genome_to_records

# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+", help="Multiple input JSONL bubble datasets.")
    ap.add_argument("--start", type=float, default=1.0)
    ap.add_argument("--stop", type=float, default=21.0)
    ap.add_argument("--step", type=float, default=0.5)
    ap.add_argument("--cap", type=float, default=21.0)
    args = ap.parse_args()

    # expand per-record decisions
    all_decisions = []
    file_mapping = []  # map each decision to file path
    for path in args.files:
        expanded = load_expand_decisions([path])
        all_decisions.extend(expanded)
        file_mapping.extend([path] * len(expanded))

    print(f"Loaded {len(all_decisions)} decision records from {len(args.files)} files.")

    # group decisions by genome
    genome_groups = defaultdict(list)
    for rec, p in zip(all_decisions, file_mapping):
        gid = genome_id_from_path(p)
        genome_groups[gid].append(rec)

    thresholds = frange(args.start, args.stop, args.step)

    # -------------------------------------------------------
    # Per-genome reporting
    # -------------------------------------------------------
    print("\n=== PER-GENOME RESULTS ===")
    for gid, recs in genome_groups.items():
        print(f"\nGenome: {gid}")
        # threshold curve
        results = evaluate_thresholds(recs, thresholds, args.cap)
        for r in results:
            acc = "nan" if math.isnan(r["accuracy"]) else f"{100*r['accuracy']:.2f}%"
            print(f"  tau={r['threshold']:>4.1f} | n={r['eval_count']:5d} | acc={acc}")

        # bubble summary
        bw, be, bc, ba = bubble_level(recs)
        if math.isnan(ba):
            print(f"  bubble_acc = nan (no evaluable bubbles)")
        else:
            print(f"  bubble_acc = {ba:.3f} ({bc}/{be})")

    # -------------------------------------------------------
    # Global report
    # -------------------------------------------------------
    print("\n=== GLOBAL SUMMARY ===")
    global_results = evaluate_thresholds(all_decisions, thresholds, args.cap)
    for r in global_results:
        acc = "nan" if math.isnan(r["accuracy"]) else f"{100*r['accuracy']:.2f}%"
        print(f"  tau={r['threshold']:>4.1f} | n={r['eval_count']:5d} | acc={acc}")

    bw, be, bc, ba = bubble_level(all_decisions)
    if math.isnan(ba):
        print("global bubble acc = nan")
    else:
        print(f"global bubble acc = {ba:.3f} ({bc}/{be})")


if __name__ == "__main__":
    main()