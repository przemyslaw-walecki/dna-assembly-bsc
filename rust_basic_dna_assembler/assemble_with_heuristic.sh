./target/release/rust_basic_dna_assembler \
  -1 ../genomes/ecoli_clean_cov20_100_300_1.fq \
  -2 ../genomes/ecoli_clean_cov20_100_300_2.fq \
  -k 21 \
  -o heurstic_assembled_ecoli_k21_cov20.fa \
  --bubble-resolver coverage \
  --bubblegun-bin BubbleGun \
  --gfa-out work/ecoli_21_graph.gfa \
  --skip-gnn \
  --bubbles-jsonl work/bubbles.json