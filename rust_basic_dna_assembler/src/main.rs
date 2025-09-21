//! Command-line interface for the AI DNA assembler.
//!
//! Loads reads from paired FASTQ files, constructs a k-mer de Bruijn graph,
//! (optionally) cleans the graph, assembles maximal unitigs, and writes:
//! - FASTA contigs to `data/{output}`
//! - GFA with KC/EC to `{output}.gfa`

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

fn main() -> io::Result<()> {
    let args = Args::parse();

    // Load reads from both FASTQ files
    let mut reads = FastqParser::new(&args.reads1).load_reads()?;
    reads.extend(FastqParser::new(&args.reads2).load_reads()?);

    // Build the k-mer graph
    let mut graph = KmerGraph::new(args.kmer);
    graph.build(&reads);

    // Optional cleaning (kept as dead code; enable as needed)
    // graph.filter_low_coverage(args.threshold);
    // graph.remove_dead_ends(Some(args.depth));
    // graph.remove_bubbles(args.bubble_depth);

    // Write GFA snapshot (KC/EC) beside FASTA output
    // let gfa_path = format!("{}.gfa", args.output);
    // graph.write_gfa(&gfa_path).unwrap();

    // Basic graph stats
    eprintln!("Loaded {} reads", reads.len());
    let edge_cnt = graph.edge_count();
    let node_cnt = graph.in_degree.len();
    let zero_indeg = graph
        .in_degree
        .iter()
        .filter(|&(_, &d)| d == 0)
        .count();
    eprintln!(
        "Graph: {} edges, {} nodes; {} sources (in_deg=0)",
        edge_cnt, node_cnt, zero_indeg
    );


    let contigs = Assembler::new(&graph).assemble_unitigs();

    // Write contigs to output FASTA in data/
    fs::create_dir_all("data")?;
    let out_path = Path::new("data").join(&args.output);
    let mut f = fs::File::create(&out_path)?;
    for (i, c) in contigs.iter().enumerate() {
        writeln!(f, ">unitig_{}", i + 1)?;
        writeln!(f, "{}", c)?;
    }
    println!("Wrote {} unitigs to {:?}", contigs.len(), out_path);

    Ok(())
}
