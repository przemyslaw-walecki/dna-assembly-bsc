 #!/bin/sh

 RUST_BACKTRACE=1 ./target/release/rust_basic_dna_assembler -1 ../genomes/klebsiella_clean_100_300_1.fq -2 ../genomes/klebsiella_clean_100_300_2.fq -k 21 -t 0 -d 0 -b 0 -o klebsiella_21_graph.gfa
 