#!/bin/bash

set -e

# Build assembler
cargo build --release

GFA_DIR="../rust_basic_dna_assembler/ecoli_gfa_out/ecoli_added"
REF_DIR="../genomes/ecoli_variants/ecoli_variants_added"
OUT_DIR="./datasets/ecoli_variants/ecoli_added"

mkdir -p "$OUT_DIR"

echo "[*] Szukam grafów w: $GFA_DIR"

# Mapowanie nazwa -> nazwa referencyjnego FASTA
declare -A REF_MAP


REF_MAP["ecoli_variant4"]="Ecoli_variant4.fasta"
REF_MAP["ecoli_variant5"]="Ecoli_variant5.fasta"

echo "[*] Rozpoczynam generowanie datasetów…"

for gfa in "$GFA_DIR"/*.gfa; do
    base=$(basename "$gfa" .gfa)

    # pierwszy człon przed _
    species="${base#*_*_}"

    echo ""
    echo "[*] Przetwarzam: $gfa (gatunek = $species)"

    # bubbles JSON:
    bubbles_json="${GFA_DIR}/${species}_bubbles.json"
    if [ ! -f "$bubbles_json" ]; then
        echo "[WARN] Brak pliku bubble JSON: $bubbles_json — pomijam"
        continue
    fi

    # odpowiedni plik referencyjny
    ref=${REF_MAP[$species]}
    if [ -z "$ref" ]; then
        echo "[ERROR] Brak mapowania referencyjnego FASTA dla: $species"
        echo "Dodaj species do REF_MAP."
        continue
    fi

    ref_path="$REF_DIR/$ref"
    if [ ! -f "$ref_path" ]; then
        echo "[ERROR] Brak referencyjnego FASTA: $ref_path"
        continue
    fi

    # output dataset
    out_jsonl="${OUT_DIR}/${species}_dataset.jsonl"

    echo "[*] bubble_json  = $bubbles_json"
    echo "[*] ref_fasta    = $ref_path"
    echo "[*] output       = $out_jsonl"

    # run dataset builder
    RUST_MIN_STACK=33554432 RAYON_NUM_THREADS=4 \
     ./target/release/rust_bubble_labelling \
        "$gfa" "$bubbles_json" "$ref_path" "$out_jsonl"

done

echo ""
echo "[DONE] Datasety wygenerowane."
