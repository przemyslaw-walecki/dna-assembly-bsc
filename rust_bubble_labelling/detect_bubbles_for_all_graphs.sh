#!/bin/sh

GFA_DIR="../rust_basic_dna_assembler/ecoli_gfa_out/ecoli_added"

echo "[*] Szukanie plików .gfa w $GFA_DIR"

for gfa in "$GFA_DIR"/*.gfa; do
    # Pomijamy gdy nie ma żadnych plików
    [ -e "$gfa" ] || { echo "Brak plików .gfa"; exit 1; }

    base=$(basename "$gfa" .gfa)
    out_json="${GFA_DIR}/${base}_bubbles.json"

    echo "[*] BubbleGun: $gfa → $out_json"

    BubbleGun -g "$gfa" bchains --bubble_json "$out_json"

    echo "Done."
done

echo "[*] Wszystkie grafy przetworzone."
