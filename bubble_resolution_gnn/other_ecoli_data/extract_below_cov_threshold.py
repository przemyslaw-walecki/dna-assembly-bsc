import argparse
import json
import math
import os
from typing import List, Dict, Any, Optional, Iterable

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
    """
    if not isinstance(paths, list) or len(paths) != 2:
        return None
    a = average_cov_for_path(edges, paths[0])
    b = average_cov_for_path(edges, paths[1])
    if a is None or b is None or a <= 0 or b <= 0:
        return None
    big, small = (a, b) if a >= b else (b, a)
    return big / small

def iter_jsonl(path: str) -> Iterable[Dict[str, Any]]:
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            yield json.loads(line)

def write_jsonl(path: str, records: Iterable[Dict[str, Any]]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")

def filter_file(input_path: str, ratio_threshold: float) -> str:
    """
    Keep all bubbles except those with label_path == null.
    Add cov_ratio if computable, then filter to cov_ratio < threshold.
    If cov_ratio cannot be computed, the record is kept (but without filtering).
    """
    base, ext = os.path.splitext(input_path)
    output_path = f"{base}_ratio_lt_{ratio_threshold:g}.jsonl"

    kept = []

    for rec in iter_jsonl(input_path):
        if rec.get("label_path") is None:
            continue

        ratio = coverage_ratio_two_paths(rec.get("edges", []), rec.get("paths", []))
        if ratio is not None:
            rec = dict(rec)
            rec["cov_ratio"] = ratio
            if ratio < ratio_threshold:
                kept.append(rec)
        else:
            kept.append(rec)

    write_jsonl(output_path, kept)
    return output_path

def main():
    p = argparse.ArgumentParser(
        description="Keep hyperbubble records with label_path not null and coverage ratio < threshold (if computable)."
    )
    p.add_argument(
        "--threshold",
        type=float,
        required=True,
        help="Coverage ratio cutoff (e.g., 3, 5, 10, 21). Keeps bubbles with ratio < threshold. If ratio missing, record is kept."
    )
    p.add_argument(
        "files",
        nargs="+",
        help="JSONL files to process (each line is one bubble record)."
    )
    args = p.parse_args()

    for ip in args.files:
        out_path = filter_file(ip, args.threshold)
        print(f"[{ip}] -> [{out_path}]")

if __name__ == "__main__":
    main()
