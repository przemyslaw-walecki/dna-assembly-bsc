#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="."
SCRIPT="augment_bubblef_dataset.py"

# Augmentation parameters
SEQ_MUT_P=0.7
SEQ_MEAN=0.35
SEQ_STD=0.15
SEQ_MINNODES=2
SEQ_MUTS_PER_SEQ=2
TRAIN_TARGET=1000
TRAIN_RATIO=0.9
COVERAGE_JITTER=0.25

for INPUT in "$INPUT_DIR"/*_cov20_ratio_lt_3.5.jsonl; do
    [ -e "$INPUT" ] || continue

    BASENAME=$(basename "$INPUT")
    NAME="${BASENAME%.jsonl}"

    TRAIN_OUT="${INPUT_DIR}/${NAME}.train.jsonl"
    TEST_OUT="${INPUT_DIR}/${NAME}.test.jsonl"

    echo "[INFO] Processing ${INPUT}"

    # Count lines
    N=$(wc -l < "$INPUT")

    echo "[INFO]   N = $N"
    echo "[INFO]   Target train = $TRAIN_TARGET"

    python3 "$SCRIPT" \
        --input "$INPUT" \
        --train-output "$TRAIN_OUT" \
        --test-output "$TEST_OUT" \
        --train-target-count "$TRAIN_TARGET" \
        --train-ratio "$TRAIN_RATIO" \
        --include-original-train \
        --seq_mutation_prob "$SEQ_MUT_P" \
        --seq-mutate-pct-mean "$SEQ_MEAN" \
        --seq-mutate-pct-std "$SEQ_STD" \
        --seq-mutate-min-nodes "$SEQ_MINNODES" \
        --seq-mutations-per-seq "$SEQ_MUTS_PER_SEQ" \
        --coverage-jitter-frac "$COVERAGE_JITTER"

done
