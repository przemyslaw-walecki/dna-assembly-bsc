#!/usr/bin/env python3
"""
Evaluate the rule: "choose the path with higher coverage" on a JSONL dataset.

Coverage for a path = SUM over edges in that path of edge["cov_mean"]
(falling back to edge["cov_min"] if "cov_mean" is missing).

Accuracy is computed against record["label_path"] ∈ {0,1}.
Ties (equal coverage) are counted separately and excluded from accuracy.
"""

import json
import math
import argparse
from statistics import mean, median

def path_coverage_sum(rec, cov_key="cov_mean", fallback_key="cov_min"):
    """Return [coverage_path0, coverage_path1] by summing edge coverages."""
    edges = rec.get("edges", [])
    paths = rec.get("paths", [])
    covs = []
    for edge_idxs in paths:
        s = 0.0
        for ei in edge_idxs:
            if not (isinstance(ei, int) and 0 <= ei < len(edges)):
                continue
            e = edges[ei]
            v = e.get(cov_key, None)
            if v is None and fallback_key is not None:
                v = e.get(fallback_key, 0)
            try:
                s += float(v)
            except Exception:
                # If it wasn't numeric, ignore that edge's coverage.
                pass
        covs.append(s)
    return covs

def evaluate_jsonl(path, cov_key="cov_mean", fallback_key="cov_min"):
    total = 0
    correct = 0
    wrong = 0
    ties = 0
    missing = 0
    deltas = []  # abs gap for non-ties

    with open(path, "r", encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                # Skip bad JSON lines
                missing += 1
                continue

            # Basic schema checks
            if "paths" not in rec or "edges" not in rec:
                missing += 1
                continue
            paths = rec["paths"]
            if not isinstance(paths, list) or len(paths) != 2:
                missing += 1
                continue
            label = rec.get("label_path", None)
            if not isinstance(label, int) or label not in (0, 1):
                missing += 1
                continue

            covs = path_coverage_sum(rec, cov_key=cov_key, fallback_key=fallback_key)
            if len(covs) != 2 or not all(math.isfinite(c) for c in covs):
                missing += 1
                continue

            # Decide by higher coverage
            if covs[0] == covs[1]:
                ties += 1
                continue

            pred = 0 if covs[0] > covs[1] else 1
            total += 1
            if pred == label:
                correct += 1
            else:
                wrong += 1
            deltas.append(abs(covs[1] - covs[0]))

    acc = (correct / total) if total else float("nan")
    return {
        "total_evaluated": total,
        "correct": correct,
        "wrong": wrong,
        "ties": ties,
        "missing_or_unusable": missing,
        "accuracy_excl_ties": acc
    }

def main():
    ap = argparse.ArgumentParser(description="Check higher-coverage rule on JSONL bubbles.")
    ap.add_argument("jsonl", help="Path to dataset (JSON Lines).")
    ap.add_argument("--cov-key", default="cov_mean",
                    help="Edge coverage field to use (default: cov_mean).")
    ap.add_argument("--fallback-key", default="cov_min",
                    help="Fallback edge coverage field if cov-key missing (default: cov_min). "
                         "Use '' to disable fallback.")
    args = ap.parse_args()

    fallback = args.fallback_key if args.fallback_key != "" else None
    stats = evaluate_jsonl(args.jsonl, cov_key=args.cov_key, fallback_key=fallback)

    # Pretty print
    def fmt(x):
        if isinstance(x, float):
            return f"{x:.6g}"
        return str(x)

    print("=== Higher-Coverage Rule Evaluation ===")
    for k in [
        "total_evaluated",
        "correct",
        "wrong",
        "ties",
        "missing_or_unusable",
        "accuracy_excl_ties"
    ]:
        print(f"{k}: {fmt(stats[k])}")

if __name__ == "__main__":
    main()
