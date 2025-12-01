#!/bin/bash

set -e

# Build assembler
cargo build --release

GFA_DIR="../rust_basic_dna_assembler/gfa_out"
REF_DIR="../genomes"
OUT_DIR="./datasets"

mkdir -p "$OUT_DIR"

echo "[*] Szukam grafów w: $GFA_DIR"

# Mapowanie nazwa -> nazwa referencyjnego FASTA
declare -A REF_MAP

# REF_MAP["Bacillus"]="Bacillus-subtilis.fasta"
# REF_MAP["Citrobacter"]="Citrobacter-freundii.fasta"
#REF_MAP["Cronobacter"]="Cronobacter-sakazakii.fasta"
#REF_MAP["Haemophilus"]="Haemophilus-influenzae.fasta"
#REF_MAP["Morganella"]="Morganella-morganii.fasta"
#REF_MAP["Proteus"]="Proteus-mirabilis.fasta"
REF_MAP["Serratia"]="Serratia-marcescens.fasta"

echo "[*] Rozpoczynam generowanie datasetów…"

for gfa in "$GFA_DIR"/*.gfa; do
    base=$(basename "$gfa" .gfa)

    # pierwszy człon przed _
    species=$(echo "$base" | cut -d'_' -f1)

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
     ./target/release/rust_basic_dna_assembler \
        "$gfa" "$bubbles_json" "$ref_path" "$out_jsonl"

done

echo ""
echo "[DONE] Datasety wygenerowane."
