
cargo clean
cargo build --release

set RAYON_NUM_THREADS=4
set RUST_MIN_STACK=33554432
.\target\release\rust_bubble_labelling.exe proteus_k21_cov10_graph.gfa proteus_k21_cov10_bubbles.json ..\genomes\Proteus-mirabilis.fasta .\datasets\proteus_dataset_k21_cov10.jsonl
.\target\release\rust_bubble_labelling.exe proteus_k21_cov30_graph.gfa proteus_k21_cov30_bubbles.json ..\genomes\Proteus-mirabilis.fasta .\datasets\proteus_dataset_k21_cov30.jsonl