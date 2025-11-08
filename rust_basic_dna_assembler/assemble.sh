 #!/bin/sh
cargo clean 
cargo build --release
 ./target/release/rust_basic_dna_assembler -1 ../genomes/klebsiella_clean_cov10_100_300_1.fq -2 ../genomes/klebsiella_clean_cov10_100_300_2.fq -k 21 -t 0 -d 0 -o klebsiella_k21_cov10_graph.gfa
 