#!/usr/bin/env python3
"""
Split-then-augment for bubblef De Bruijn GNN datasets.

Flow:
1) Load originals from --input (JSONL).
2) Deterministically split originals into train/test by --train-ratio.
3) Write test set *untouched* to --test-output.
4) Optionally write original train set (use --include-original-train).
5) Augment ONLY the train set until --train-target-count is reached.
   - Sequence augmentation: mutate a Gaussian-sampled percentage of sequences
     (node sequences, and optionally start_seq/end_seq) with a fixed number of base flips per sequence.
   - Coverage augmentation: multiplicative jitter of node/edge coverage values.
   - New items get bubble_id reassigned and metadata: augmented, aug_parent, aug_method.

Notes:
- Writes JSONL for each split. The test file is always re-written from originals.
- Train file is appended in --resume mode; progress is fsynced periodically.

Augmentation knobs:
- --seq_mutation_prob (default 0.6): per-item probability to use sequence mutation; otherwise coverage jitter is used.
- Sequence mutation controls:
    --seq-mutate-pct-mean (default 0.2)
    --seq-mutate-pct-std  (default 0.1)
    --seq-mutate-min-nodes (default 1)
    --seq-include-terminal-prob (default 0.2)
    --seq-mutations-per-seq (default 1)
- Coverage jitter controls:
    --coverage-jitter-frac (default 0.05)
"""
import argparse
import json
import math
import os
import random
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

BASES = ("A", "C", "G", "T")

def load_jsonl(path: Path) -> List[Dict[str, Any]]:
    data = []
    with path.open() as f:
        for ln, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
                data.append(obj)
            except Exception as e:
                print(f"[WARN] Skipping line {ln}: JSON error: {e}", file=sys.stderr, flush=True)
    return data

def save_jsonl(path: Path, items: List[Dict[str, Any]]) -> None:
    with path.open("w") as f:
        for obj in items:
            f.write(json.dumps(obj, ensure_ascii=False) + "\n")

def read_count_and_max_id(path: Path) -> Tuple[int,int]:
    if not path.exists():
        return 0, 0
    cnt = 0
    max_id = 0
    with path.open("r") as f:
        for line in f:
            cnt += 1
            try:
                o = json.loads(line)
                if "bubble_id" in o:
                    max_id = max(max_id, int(o["bubble_id"]))
            except Exception:
                pass
    return cnt, max_id

def rand_base_except(b: str) -> str:
    b = b.upper()
    choices = [x for x in BASES if x != b]
    return random.choice(choices) if choices else random.choice(BASES)

def mutate_one_seq(seq: str, n_mut: int = 1) -> str:
    if not seq or n_mut <= 0:
        return seq
    s = list(seq)
    k = min(n_mut, len(s))
    if k <= 0:
        return seq
    positions = random.sample(range(len(s)), k=k)
    for pos in positions:
        s[pos] = rand_base_except(s[pos])
    return "".join(s)

def jitter_value(val, jitter_frac: float):
    factor = 1.0 + random.uniform(-jitter_frac, jitter_frac)
    new_val = val * factor
    if isinstance(val, int):
        return max(0, int(round(new_val)))
    if isinstance(val, float):
        return float(new_val)
    return val

def deep_copy(obj: Dict[str, Any]) -> Dict[str, Any]:
    return json.loads(json.dumps(obj))

def _clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))

def _sample_gaussian_pct(mean: float, std: float) -> float:
    # Clamp to [0, 1] after sampling
    return _clamp(random.gauss(mean, std), 0.0, 1.0)

def seq_mutation_gaussian(
    item: Dict[str, Any],
    pct_mean: float = 0.2,
    pct_std: float = 0.1,
    min_nodes: int = 1,
    include_terminal_prob: float = 0.2,
    mutations_per_seq: int = 1,
) -> Tuple[Dict[str, Any], str]:
    """
    Mutate a Gaussian-sampled percentage of sequences in the item.
    Pool = all node sequences; with probability include_terminal_prob, also include start_seq/end_seq if present.
    Each selected sequence gets 'mutations_per_seq' base flips.
    """
    out = deep_copy(item)

    # Build mutation pool: tuples of (kind, index_or_key, current_seq_string)
    # kind in {"node", "start_seq", "end_seq"}
    pool: List[Tuple[str, Any, str]] = []

    nodes = out.get("nodes", [])
    if isinstance(nodes, list):
        for i, n in enumerate(nodes):
            if isinstance(n, dict) and isinstance(n.get("seq"), str) and n["seq"]:
                pool.append(("node", i, n["seq"]))

    # Optionally include terminals as pseudo-nodes
    if random.random() < include_terminal_prob:
        if isinstance(out.get("start_seq"), str) and out["start_seq"]:
            pool.append(("start_seq", "start_seq", out["start_seq"]))
        if isinstance(out.get("end_seq"), str) and out["end_seq"]:
            pool.append(("end_seq", "end_seq", out["end_seq"]))

    if not pool:
        return out, "no_op_seq_mutation"

    # Sample percentage and compute how many unique sequences to mutate
    pct = _sample_gaussian_pct(pct_mean, pct_std)
    n_pool = len(pool)
    n_to_mutate = max(0, int(round(pct * n_pool)))
    if n_to_mutate <= 0 and min_nodes > 0:
        n_to_mutate = min(min_nodes, n_pool)

    if n_to_mutate <= 0:
        return out, "no_op_seq_mutation"

    # Choose disjoint items from pool
    chosen = random.sample(pool, k=n_to_mutate)

    # Map of old_seq -> new_seq (for edge consistency updates)
    # Note: if the same exact string appears multiple times, they will share a single mapping.
    seq_map: Dict[str, str] = {}

    mutated_count = 0
    for kind, idx_or_key, old_seq in chosen:
        new_seq = mutate_one_seq(old_seq, n_mut=mutations_per_seq)
        if new_seq == old_seq:
            # could happen if seq empty or n_mut==0
            continue
        seq_map[old_seq] = new_seq
        if kind == "node":
            out["nodes"][idx_or_key]["seq"] = new_seq
        elif kind in ("start_seq", "end_seq"):
            out[idx_or_key] = new_seq
        mutated_count += 1

    # Propagate sequence renames to edges
    if mutated_count > 0 and isinstance(out.get("edges"), list):
        for e in out["edges"]:
            if isinstance(e, dict):
                if "source_seq" in e and isinstance(e["source_seq"], str):
                    if e["source_seq"] in seq_map:
                        e["source_seq"] = seq_map[e["source_seq"]]
                if "target_seq" in e and isinstance(e["target_seq"], str):
                    if e["target_seq"] in seq_map:
                        e["target_seq"] = seq_map[e["target_seq"]]

    if mutated_count == 0:
        return out, "no_op_seq_mutation"

    # Encode method with details for traceability
    method = f"seq_mutation_gaussian:pct={pct:.3f},mutated={mutated_count}/{n_pool},muts_per_seq={mutations_per_seq}"
    return out, method

def coverage_jitter(item: Dict[str, Any], jitter_frac: float = 0.05) -> Tuple[Dict[str, Any], str]:
    out = deep_copy(item)
    method = "coverage_jitter"

    if isinstance(out.get("nodes"), list):
        for n in out["nodes"]:
            if isinstance(n, dict) and "cov" in n and isinstance(n["cov"], (int, float)):
                n["cov"] = jitter_value(n["cov"], jitter_frac)

    if isinstance(out.get("edges"), list):
        for e in out["edges"]:
            if isinstance(e, dict):
                if "cov_mean" in e and isinstance(e["cov_mean"], (int, float)):
                    e["cov_mean"] = jitter_value(e["cov_mean"], jitter_frac)
                if "cov_min" in e and isinstance(e["cov_min"], (int, float)):
                    e["cov_min"] = jitter_value(e["cov_min"], jitter_frac)
                if ("cov_mean" in e and "cov_min" in e
                    and isinstance(e["cov_mean"], (int,float))
                    and isinstance(e["cov_min"], (int,float))):
                    if e["cov_min"] > e["cov_mean"]:
                        e["cov_min"] = max(0, int(e["cov_mean"]))
    return out, method

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--train-output", required=True, type=Path)
    ap.add_argument("--test-output", required=True, type=Path)
    ap.add_argument("--train-target-count", type=int, required=True,
                    help="Total train output lines (originals + augmentations if --include-original-train).")
    ap.add_argument("--train-ratio", type=float, default=0.8,
                    help="Fraction of originals to use for train before augmentation.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--include-original-train", action="store_true",
                    help="Write original train items before augmentations.")
    # Sequence mutation knobs
    ap.add_argument("--seq_mutation_prob", type=float, default=0.6,
                    help="Per-item probability to apply sequence mutation instead of coverage jitter.")
    ap.add_argument("--seq-mutate-pct-mean", type=float, default=0.2,
                    help="Mean fraction of sequences to mutate (clamped to [0,1]).")
    ap.add_argument("--seq-mutate-pct-std", type=float, default=0.1,
                    help="Stddev of fraction of sequences to mutate.")
    ap.add_argument("--seq-mutate-min-nodes", type=int, default=1,
                    help="Minimum number of sequences to mutate when pool not empty.")
    ap.add_argument("--seq-include-terminal-prob", type=float, default=0.2,
                    help="Probability to include start_seq/end_seq in mutation pool.")
    ap.add_argument("--seq-mutations-per-seq", type=int, default=1,
                    help="Number of base flips per selected sequence.")
    # Coverage jitter
    ap.add_argument("--coverage-jitter-frac", type=float, default=0.05)
    # Misc
    ap.add_argument("--progress-interval", type=int, default=10000)
    ap.add_argument("--resume", action="store_true",
                    help="Append to existing train-output until train-target-count is reached. Test output is always rewritten from originals.")
    args = ap.parse_args()

    random.seed(args.seed)

    base = load_jsonl(args.input)
    if not base:
        print("[ERROR] No valid input.", file=sys.stderr, flush=True)
        sys.exit(2)

    # deterministic shuffle before split (seeded)
    rng = random.Random(args.seed)
    idxs = list(range(len(base)))
    rng.shuffle(idxs)
    split_idx = int(len(base) * args.train_ratio)
    train_idx = idxs[:split_idx]
    test_idx  = idxs[split_idx:]

    train_base = [base[i] for i in train_idx]
    test_base  = [base[i] for i in test_idx]

    # 1) Always write test set clean (no augmentation)
    save_jsonl(args.test_output, test_base)
    print(f"[INFO] Wrote test originals: {len(test_base)} -> {args.test_output}", file=sys.stderr, flush=True)

    # 2) Prepare train output (resume or fresh)
    if args.resume:
        existing, max_out_id = read_count_and_max_id(args.train_output)
        mode = "a" if args.train_output.exists() else "w"
        out_f = args.train_output.open(mode)
        written = existing
        print(f"[INFO] Resume train. Existing: {existing}, max bubble_id in train file: {max_out_id}", file=sys.stderr, flush=True)
    else:
        out_f = args.train_output.open("w")
        written = 0
        max_out_id = 0

    # include original train samples (once) at the top when not resuming or empty file
    if args.include_original_train and written == 0:
        for obj in train_base:
            out_f.write(json.dumps(obj, ensure_ascii=False) + "\n")
            written += 1
        out_f.flush()
        os.fsync(out_f.fileno())
        print(f"[INFO] Wrote train originals: {len(train_base)}", file=sys.stderr, flush=True)

    # figure next bubble_id to assign for augmented items
    try:
        max_orig_id = max(o.get("bubble_id", 0) for o in base)
    except Exception:
        max_orig_id = 0
    next_id = max(max_out_id + 1, max_orig_id + 1)

    # 3) Augment only from train_base until reaching train-target-count
    i = 0
    while written < args.train_target_count:
        src = train_base[i % len(train_base)]
        parent_id = src.get("bubble_id")

        if random.random() < args.seq_mutation_prob:
            aug, method = seq_mutation_gaussian(
                src,
                pct_mean=args.seq_mutate_pct_mean,
                pct_std=args.seq_mutate_pct_std,
                min_nodes=args.seq_mutate_min_nodes,
                include_terminal_prob=args.seq_include_terminal_prob,
                mutations_per_seq=args.seq_mutations_per_seq,
            )
        else:
            aug, method = coverage_jitter(src, jitter_frac=args.coverage_jitter_frac)

        # assign metadata
        aug["bubble_id"] = next_id
        aug["augmented"] = True
        aug["aug_parent"] = parent_id
        aug["aug_method"] = method

        out_f.write(json.dumps(aug, ensure_ascii=False) + "\n")
        written += 1
        next_id += 1
        i += 1

        if written % args.progress_interval == 0 or written == args.train_target_count:
            out_f.flush()
            os.fsync(out_f.fileno())
            print(f"[INFO] Generated train {written} / {args.train_target_count}...", file=sys.stderr, flush=True)

    out_f.close()
    print(f"[INFO] Done. Train lines: {written}. Test lines: {len(test_base)}.", file=sys.stderr, flush=True)

if __name__ == "__main__":
    main()
