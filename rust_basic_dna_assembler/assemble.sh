 #!/bin/sh
cargo clean 
cargo build --release
./target/release/rust_basic_dna_assembler -1 ../genomes/ecoli_clean_cov10_100_300_1.fq -2 ../genomes/ecoli_clean_cov10_100_300_2.fq -k 21 -t 0 -d 0 -o ecoli_k21_cov10_graph.gfa
./target/release/rust_basic_dna_assembler -1 ../genomes/ecoli_clean_cov30_100_300_1.fq -2 ../genomes/ecoli_clean_cov30_100_300_2.fq -k 21 -t 0 -d 0 -o ecoli_k21_cov30_graph.gfa
./target/release/rust_basic_dna_assembler -1 ../genomes/klebsiella_clean_cov30_100_300_1.fq -2 ../genomes/klebsiella_clean_cov30_100_300_2.fq -k 21 -t 0 -d 0 -o klebsiella_k21_cov30_graph.gfa
./target/release/rust_basic_dna_assembler -1 ../genomes/salmonella_clean_cov30_100_300_1.fq -2 ../genomes/salmonella_clean_cov30_100_300_2.fq -k 21 -t 0 -d 0 -o salmonella_k21_cov30_graph.gfa