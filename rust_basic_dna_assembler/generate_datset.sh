 #!/bin/sh

cargo clean
cargo build --release

RAYON_NUM_THREADS=4 ./target/release/rust_basic_dna_assembler ecoli_21_graph.gfa ecoli_k21_bubbles.json ../genomes/Escherichia_coli_str_K-12_substr_MG1655.fasta ./ecoli_dataset_new.jsonl 
RAYON_NUM_THREADS=4 ./target/release/rust_basic_dna_assembler klebsiella_21_graph.gfa klebsiella_k21_bubbles.json ../genomes/Klebsiella_pneumoniae_subsp_pneumoniae_NTUH-K2044_DNA.fasta ./klebsiella_dataset_new.jsonl