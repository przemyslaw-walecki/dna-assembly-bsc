# Salmonella

./target/release/rust_basic_dna_assembler \
  -1 ../genomes/proteus_clean_cov20_100_300_1.fq \
  -2 ../genomes/proteus_clean_cov20_100_300_2.fq \
  -k 21 \
  -o ai_assembled_proteus_k21_cov20_04_01.fa \
  --gfa-out work/proteus_21_graph.gfa \
  --bubblegun-bin BubbleGun \
  --bubbles-jsonl work/proteus_k21_bubbles.json \
  --python-bin /usr/bin/python3 \
  --gnn-script ../bubble_resolution_gnn/GNN_infer.py \
  --decisions-jsonl work/proteus_k21_decisions.jsonl \
  --gnn-args='--gfa work/proteus_21_graph.gfa --ckpt ../bubble_resolution_gnn/architecture_sweep_24_12/models/emb8_g16_e256_mean_cov20_fold3_Proteus_run3.pth --cpu'

./target/release/rust_basic_dna_assembler \
  -1 ../genomes/proteus_clean_cov20_100_300_1.fq \
  -2 ../genomes/proteus_clean_cov20_100_300_2.fq \
  -k 21 \
  -o heurstic_assembled_proteus_k21_cov20_04_01.fa \
  --bubble-resolver coverage \
  --bubblegun-bin BubbleGun \
  --gfa-out work/proteus_21_graph.gfa \
  --skip-gnn \
  --bubbles-jsonl work/bubbles.json

# Klebsiella 

./target/release/rust_basic_dna_assembler \
  -1 ../genomes/klebsiella_clean_cov20_100_300_1.fq \
  -2 ../genomes/klebsiella_clean_cov20_100_300_2.fq \
  -k 21 \
  -o ai_assembled_klebsiella_k21_cov20_04_01.fa \
  --gfa-out work/klebsiella_21_graph.gfa \
  --bubblegun-bin BubbleGun \
  --bubbles-jsonl work/klebsiella_k21_bubbles.json \
  --python-bin /usr/bin/python3 \
  --gnn-script ../bubble_resolution_gnn/GNN_infer.py \
  --decisions-jsonl work/klebsiella_k21_decisions.jsonl \
  --gnn-args='--gfa work/klebsiella_21_graph.gfa --ckpt ../bubble_resolution_gnn/architecture_sweep_24_12/models/emb8_g16_e256_mean_cov20_fold2_Klebsiella_run3.pth --cpu'

./target/release/rust_basic_dna_assembler \
  -1 ../genomes/klebsiella_clean_cov20_100_300_1.fq \
  -2 ../genomes/klebsiella_clean_cov20_100_300_2.fq \
  -k 21 \
  -o heurstic_assembled_klebsiella_k21_cov20_04_01.fa \
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
  -o ai_assembled_ecoli_k21_cov20_04_01.fa \
  --gfa-out work/ecoli_21_graph.gfa \
  --bubblegun-bin BubbleGun \
  --bubbles-jsonl work/ecoli_k21_bubbles.json \
  --python-bin /usr/bin/python3 \
  --gnn-script ../bubble_resolution_gnn/GNN_infer.py \
  --decisions-jsonl work/ecoli_k21_decisions.jsonl \
  --gnn-args='--gfa work/ecoli_21_graph.gfa --ckpt ../bubble_resolution_gnn/architecture_sweep_24_12/models/emb8_g16_e256_mean_cov20_fold1_Ecoli_run3.pth --cpu'

./target/release/rust_basic_dna_assembler \
  -1 ../genomes/ecoli_clean_cov20_100_300_1.fq \
  -2 ../genomes/ecoli_clean_cov20_100_300_2.fq \
  -k 21 \
  -o heurstic_assembled_ecoli_k21_cov20_04_01.fa \
  --bubble-resolver coverage \
  --bubblegun-bin BubbleGun \
  --gfa-out work/ecoli_21_graph.gfa \
  --skip-gnn \
  --bubbles-jsonl work/bubbles.json


