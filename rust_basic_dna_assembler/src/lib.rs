//! Biblioteczna warstwa crate: wystawia moduły z `src/` do użycia w testach
//! integracyjnych (`tests/*.rs`) oraz (opcjonalnie) przez inne craty.

pub mod fastq_parser;
pub mod kmer_graph;
pub mod assembler;

// Wygodne re-exporty (opcjonalne, ale ułatwiają import w testach).
pub use assembler::Assembler;
pub use fastq_parser::FastqParser;
pub use kmer_graph::KmerGraph;
