
cargo clean
cargo build --release

set RAYON_NUM_THREADS=4
set RUST_MIN_STACK=33554432
.\target\release\rust_basic_dna_assembler.exe klebsiella_k21_cov30_graph.gfa klebsiella_k21_cov30_bubbles.json ..\genomes\Klebsiella_pneumoniae_subsp_pneumoniae_NTUH-K2044_DNA.fasta .\klebsiella_dataset_k21_cov30.jsonl