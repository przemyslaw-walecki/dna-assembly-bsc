#!/usr/bin/env bash
set -euo pipefail

ASSEMBLER="./target/release/rust_basic_dna_assembler"
READ_DIR="../genomes/ecoli_variants/ecoli_variants_added"
OUT_DIR="./ecoli_gfa_out"
K=21

mkdir -p "$OUT_DIR"

for f1 in "$READ_DIR"/*_1.fq; do
    base=$(basename "$f1" "_clean_cov20_100_300_1.fq")
    f2="$READ_DIR/${base}_clean_cov20_100_300_2.fq"

    if [[ ! -f "$f2" ]]; then
        echo "[warn] Missing pair for $f1, skipping"
        continue
    fi

    gfa_out="$OUT_DIR/${base}.gfa"

    echo ">>> Building graph for $base"

    $ASSEMBLER \
        -1 "$f1" \
        -2 "$f2" \
        -k $K \
        --skip-bubblegun \
        --skip-gnn \
        --bubble-resolver none \
        --gfa-out "$gfa_out" \
        --output /dev/null
done
