#!/usr/bin/env python3
"""
Split-then-augment for bubblef De Bruijn GNN datasets.

Flow:
1) Load originals from --input (JSONL).
2) Deterministically split originals into train/test by --train-ratio.
3) Write test set *untouched* to --test-output.
4) Optionally write original train set (use --include-original-train).
5) Augment ONLY the train set until --train-target-count is reached.
   - Light augmentations only (sequence mutation or coverage jitter).
   - New items get bubble_id reassigned and metadata: augmented, aug_parent, aug_method.

This avoids data leakage by ensuring the test split contains no augmented
variants of training samples.

Examples:
    python augment_bubblef_split_then_augment.py \
        --input ecoli_dataset_cov20_full.jsonl \
        --train-output ecoli_train_aug.jsonl \
        --test-output ecoli_test.jsonl \
        --train-target-count 320000 \
        --train-ratio 0.8 \
        --seed 42 \
        --include-original-train

    # Resume appending train augmentations if a previous run exists:
    python augment_bubblef_split_then_augment.py \
        --input ecoli_dataset_cov20_full.jsonl \
        --train-output ecoli_train_aug.jsonl \
        --test-output ecoli_test.jsonl \
        --train-target-count 320000 \
        --train-ratio 0.8 \
        --seed 42 \
        --include-original-train \
        --resume

Notes:
- Writes JSONL for each split. The test file is always re-written from originals.
- Train file is appended in --resume mode; progress is fsynced periodically.
- Augmentation knobs:
    --seq-mutations-per-item (default 1), --coverage-jitter-frac (default 0.05),
    --seq_mutation_prob (default 0.6).
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
    positions = random.sample(range(len(s)), k=min(n_mut, len(s)))
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

def seq_mutation(item: Dict[str, Any], seq_mutations_per_item: int = 1) -> Tuple[Dict[str, Any], str]:
    out = deep_copy(item)
    method = "seq_mutation"
    nodes = out.get("nodes", [])
    mutated_any = False

    if random.random() < 0.2 and ("start_seq" in out or "end_seq" in out):
        keys = [k for k in ("start_seq", "end_seq") if k in out]
        if keys:
            key = random.choice(keys)
            old_seq = out[key]
            if isinstance(old_seq, str) and old_seq:
                new_seq = mutate_one_seq(old_seq, n_mut=seq_mutations_per_item)
                out[key] = new_seq
                mutated_any = True
                for e in out.get("edges", []) or []:
                    if isinstance(e, dict):
                        if e.get("source_seq") == old_seq:
                            e["source_seq"] = new_seq
                        if e.get("target_seq") == old_seq:
                            e["target_seq"] = new_seq
    elif isinstance(nodes, list) and nodes:
        idx = random.randrange(len(nodes))
        node = nodes[idx]
        if isinstance(node, dict) and isinstance(node.get("seq"), str) and node["seq"]:
            old_seq = node["seq"]
            new_seq = mutate_one_seq(old_seq, n_mut=seq_mutations_per_item)
            node["seq"] = new_seq
            mutated_any = True
            for e in out.get("edges", []) or []:
                if isinstance(e, dict):
                    if e.get("source_seq") == old_seq:
                        e["source_seq"] = new_seq
                    if e.get("target_seq") == old_seq:
                        e["target_seq"] = new_seq

    if not mutated_any:
        method = "no_op_seq_mutation"
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
                if "cov_mean" in e and "cov_min" in e and isinstance(e["cov_mean"], (int,float)) and isinstance(e["cov_min"], (int,float)):
                    if e["cov_min"] > e["cov_mean"]:
                        e["cov_min"] = max(0, int(e["cov_mean"]))
    return out, method

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--train-output", required=True, type=Path)
    ap.add_argument("--test-output", required=True, type=Path)
    ap.add_argument("--train-target-count", type=int, required=True, help="Total train output lines (originals + augmentations if --include-original-train).")
    ap.add_argument("--train-ratio", type=float, default=0.8, help="Fraction of originals to use for train before augmentation.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--include-original-train", action="store_true", help="Write original train items before augmentations.")
    ap.add_argument("--seq-mutations-per-item", type=int, default=1)
    ap.add_argument("--coverage-jitter-frac", type=float, default=0.05)
    ap.add_argument("--seq_mutation_prob", type=float, default=0.6)
    ap.add_argument("--progress-interval", type=int, default=10000)
    ap.add_argument("--resume", action="store_true", help="Append to existing train-output until train-target-count is reached. Test output is always rewritten from originals.")
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
            aug, method = seq_mutation(src, seq_mutations_per_item=args.seq_mutations_per_item)
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
