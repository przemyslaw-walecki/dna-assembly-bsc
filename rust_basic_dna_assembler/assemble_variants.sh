# Salmonella

./target/release/rust_basic_dna_assembler \
  -1 ../genomes/salmonella_clean_cov20_100_300_1.fq \
  -2 ../genomes/salmonella_clean_cov20_100_300_2.fq \
  -k 21 \
  -o ai_assembled_salmonella_k21_cov20.fa \
  --gfa-out work/salmonella_21_graph.gfa \
  --bubblegun-bin BubbleGun \
  --bubbles-jsonl work/salmonella_k21_bubbles.json \
  --python-bin /usr/bin/python3 \
  --gnn-script ../bubble_resolution_gnn/GNN_infer.py \
  --decisions-jsonl work/salmonella_k21_decisions.jsonl \
  --gnn-args='--gfa work/salmonella_21_graph.gfa --ckpt ../bubble_resolution_gnn/cov_sweep_runs/cov20_fold3_salmonella_emb16_g64_e64_mean.pth --directml --gfa work/salmonella_21_graph.gfa' 

./target/release/rust_basic_dna_assembler \
  -1 ../genomes/salmonella_clean_cov20_100_300_1.fq \
  -2 ../genomes/salmonella_clean_cov20_100_300_2.fq \
  -k 21 \
  -o heurstic_assembled_salmonella_k21_cov20.fa \
  --bubble-resolver coverage \
  --bubblegun-bin BubbleGun \
  --gfa-out work/salmonella_21_graph.gfa \
  --skip-gnn \
  --bubbles-jsonl work/bubbles.json

# Klebsiella 

./target/release/rust_basic_dna_assembler \
  -1 ../genomes/klebsiella_clean_cov20_100_300_1.fq \
  -2 ../genomes/klebsiella_clean_cov20_100_300_2.fq \
  -k 21 \
  -o ai_assembled_klebsiella_k21_cov20.fa \
  --gfa-out work/klebsiella_21_graph.gfa \
  --bubblegun-bin BubbleGun \
  --bubbles-jsonl work/klebsiella_k21_bubbles.json \
  --python-bin /usr/bin/python3 \
  --gnn-script ../bubble_resolution_gnn/GNN_infer.py \
  --decisions-jsonl work/klebsiella_k21_decisions.jsonl \
  --gnn-args='--gfa work/klebsiella_21_graph.gfa --ckpt ../bubble_resolution_gnn/cov_sweep_runs/cov20_fold2_klebsiella_emb16_g64_e64_mean.pth --directml --gfa work/klebsiella_21_graph.gfa' 

./target/release/rust_basic_dna_assembler \
  -1 ../genomes/klebsiella_clean_cov20_100_300_1.fq \
  -2 ../genomes/klebsiella_clean_cov20_100_300_2.fq \
  -k 21 \
  -o heurstic_assembled_klebsiella_k21_cov20.fa \
  --bubble-resolver coverage \
  --bubblegun-bin BubbleGun \
  --gfa-out work/klebsiella_21_graph.gfa \
  --skip-gnn \
  --bubbles-jsonl work/bubbles.json

# Ecoli

./target/release/rust_basic_dna_assembler \
  -1 ../genomes/ecoli_clean_cov20_100_300_1.fq \
  -2 ../genomes/ecoli_clean_cov20_100_300_2.fq \
  -k 21 \
  -o ai_assembled_ecoli_k21_cov20.fa \
  --gfa-out work/ecoli_21_graph.gfa \
  --bubblegun-bin BubbleGun \
  --bubbles-jsonl work/ecoli_k21_bubbles.json \
  --python-bin /usr/bin/python3 \
  --gnn-script ../bubble_resolution_gnn/GNN_infer.py \
  --decisions-jsonl work/ecoli_k21_decisions.jsonl \
  --gnn-args='--gfa work/ecoli_21_graph.gfa --ckpt ../bubble_resolution_gnn/cov_sweep_runs/cov20_fold1_ecoli_emb16_g64_e64_mean.pth  --directml --gfa work/ecoli_21_graph.gfa' 

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


