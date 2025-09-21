 #!/bin/sh

cargo clean
cargo build --release

RAYON_NUM_THREADS=4 ./target/release/rust_basic_dna_assembler klebsiella_k21_cov20_graph.gfa klebsiella_k21_cov20_bubbles.json ../genomes/Klebsiella_pneumoniae_subsp_pneumoniae_NTUH-K2044_DNA.fasta ./klebsiella_dataset_cov20_full_upgraded.jsonl

RAYON_NUM_THREADS=4 ./target/release/rust_basic_dna_assembler salmonella_k21_cov20_graph.gfa salmonella_k21_cov20_bubbles.json ../genomes/Salmonella_enterica_subsp_enterica_serovar_Typhimurium_str_LT2.fasta ./salmonella_dataset_cov20_full_upgraded.jsonl

RAYON_NUM_THREADS=4 ./target/release/rust_basic_dna_assembler ecoli_k21_cov20_graph.gfa ecoli_k21_cov20_bubbles.json ../genomes/Escherichia_coli_str_K-12_substr_MG1655.fasta ./ecoli_dataset_cov20_full_upgraded.jsonl