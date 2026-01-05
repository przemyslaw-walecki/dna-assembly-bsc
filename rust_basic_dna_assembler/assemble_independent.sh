./target/release/rust_basic_dna_assembler \
  -1 ../genomes/salmonella_clean_cov20_100_300_1.fq \
  -2 ../genomes/salmonella_clean_cov20_100_300_2.fq \
  -k 21 \
  -o ai_assembled_salmonella_k21_cov20_04_01.fa \
  --gfa-out work/salmonella_21_graph.gfa \
  --bubblegun-bin BubbleGun \
  --bubbles-jsonl work/salmonella_k21_bubbles.json \
  --python-bin /usr/bin/python3 \
  --gnn-script ../bubble_resolution_gnn/GNN_infer.py \
  --decisions-jsonl work/salmonella_k21_decisions.jsonl \
  --gnn-args='--gfa work/salmonella_21_graph.gfa --ckpt ../bubble_resolution_gnn/train_final/gnn_no_salmonella_in_train.pth --cpu'

./target/release/rust_basic_dna_assembler \
  -1 ../genomes/salmonella_clean_cov20_100_300_1.fq \
  -2 ../genomes/salmonella_clean_cov20_100_300_2.fq \
  -k 21 \
  -o heurstic_assembled_salmonella_k21_cov20_04_01.fa \
  --bubble-resolver coverage \
  --bubblegun-bin BubbleGun \
  --gfa-out work/salmonella_21_graph.gfa \
  --skip-gnn \
  --bubbles-jsonl work/bubbles.json
