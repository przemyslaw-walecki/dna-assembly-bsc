#!/bin/sh
# cargo clean
# cargo build --release

./target/release/rust_basic_dna_assembler \
  -1 ../genomes/ecoli_clean_cov20_100_300_1.fq \
  -2 ../genomes/ecoli_clean_cov20_100_300_2.fq \
  -k 21 \
  -o "./ai_assembled_ecoli_k21_cov20.fasta" \
  --gfa-out work/ecoli_21_graph.gfa \
  --bubblegun-bin BubbleGun \
  --bubbles-jsonl work/ecoli_k21_bubbles.json \
  --python-bin /usr/bin/python3 \
  --gnn-script ../bubble_resolution_gnn/GNN_infer.py \
  --decisions-jsonl work/ecoli_k21_decisions.jsonl \
  --gnn-args='--gfa work/ecoli_21_graph.gfa --ckpt ../bubble_resolution_gnn/model2.8k.pth --directml --gfa work/ecoli_21_graph.gfa' 