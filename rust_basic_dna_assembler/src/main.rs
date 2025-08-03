//! Command-line interface for the AI DNA assembler.
//!
//! This module loads sequencing reads from paired FASTQ files, constructs a k-mer de Bruijn graph,
//! performs graph cleaning operations, assembles contigs, and writes the result to a FASTA file.
//!
//! # Functionality
//!
//! - Parses command-line arguments
//! - Loads reads using `FastqParser`
//! - Builds and cleans the graph via `KmerGraph`
//! - Assembles contigs with `Assembler`
//! - Outputs contigs to `data/{output_file}` in FASTA format

mod assembler;
mod fastq_parser;
mod kmer_graph;

use assembler::Assembler;
use clap::Parser;
use fastq_parser::FastqParser;
use kmer_graph::KmerGraph;
use std::fs;
use std::io::{self, Write};
use std::path::Path;

/// Command-line arguments for the assembler.
///
/// Provides options for input files, k-mer length, filtering threshold, pruning depth,
/// bubble detection depth, and output path.
#[derive(Parser)]
struct Args {
    /// Path to the first FASTQ file (`-1`)
    #[clap(short = '1')]
    reads1: String,

    /// Path to the second FASTQ file (`-2`)
    #[clap(short = '2')]
    reads2: String,

    /// Length of k-mers to use in the de Bruijn graph (`-k`)
    #[clap(short, long, default_value = "55")]
    kmer: usize,

    /// Minimum edge count to retain in the graph (`-t`)
    #[clap(short, long, default_value = "0")]
    threshold: usize,

    /// Maximum depth for dead-end pruning (`-d`)
    #[clap(short = 'd', long, default_value = "5")]
    depth: usize,

    /// Maximum depth for bubble removal (`-b`)
    #[clap(short = 'b', long, default_value = "20")]
    bubble_depth: usize,

    /// Output file name for contigs (`-o`)
    #[clap(short = 'o', long, default_value = "assembled.fa")]
    output: String,
}

/// Main entry point for the assembler pipeline.
///
/// This function loads the reads, constructs and cleans the graph, assembles the contigs,
/// and writes them to a FASTA file.
///
/// # Errors
/// Returns `io::Error` if reading or writing files fails.
fn main() -> io::Result<()> {
    let args = Args::parse();

    // Load reads from both FASTQ files
    let mut reads = FastqParser::new(&args.reads1).load_reads()?;
    reads.extend(FastqParser::new(&args.reads2).load_reads()?);

    // Build and clean the k-mer graph
    let mut graph = KmerGraph::new(args.kmer);

    graph.build(&reads);
    //graph.filter_low_coverage(args.threshold);
    // graph.remove_dead_ends(Some(args.depth));
    graph.write_gfa(&args.output).unwrap();
    // graph.remove_bubbles(args.bubble_depth);

    // Basic graph stats
    eprintln!("Loaded {} reads", reads.len());
    let edge_cnt = graph.edge_count();
    let node_cnt = graph.in_degree.len();
    let zero_indeg: Vec<_> = graph
        .in_degree
        .iter()
        .filter_map(|(&u, &d)| if d == 0 { Some(u) } else { None })
        .collect();
    eprintln!(
        "Graph: {} edges, {} nodes; {} sources (in_deg=0)",
        edge_cnt, node_cnt, zero_indeg.len()
    );
    //for &u in zero_indeg.iter().take(5) {
    //    eprintln!(" source kmer: {}", graph.decode(u));
    //}

    // Assemble contigs from the cleaned graph
    //let contigs = Assembler::new(&graph).assemble_contigs();

    Write contigs to output file in FASTA format
    fs::create_dir_all("data")?;
    let out_path = Path::new("data").join(&args.output);
    let mut f = fs::File::create(&out_path)?;
    for (i, c) in contigs.iter().enumerate() {
        writeln!(f, ">contig_{}", i + 1)?;
        writeln!(f, "{}", c)?;
    }

    // println!("Wrote {} contigs to {:?}", contigs.len(), out_path);
    Ok(())
}
