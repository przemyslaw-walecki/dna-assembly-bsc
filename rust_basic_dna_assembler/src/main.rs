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

#[derive(Parser)]
struct Args {
    #[clap(short = '1')]
    reads1: String,
    #[clap(short = '2')]
    reads2: String,
    #[clap(short, long, default_value = "55")]
    kmer: usize,
    #[clap(short, long, default_value = "5")]
    threshold: usize,
    #[clap(short = 'd', long, default_value = "5")]
    depth: usize,
    #[clap(short = 'o', long, default_value = "assembled.fa")]
    output: String,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    // Load reads
    let mut reads = FastqParser::new(&args.reads1).load_reads()?;
    reads.extend(FastqParser::new(&args.reads2).load_reads()?);

    // Build graph
    let mut graph = KmerGraph::new(args.kmer);
    graph.build(&reads);
    graph.filter_low_coverage(args.threshold);
    graph.remove_dead_ends(Some(args.depth));
    graph.remove_bubbles(20);

    // Assemble contigs
    let contigs = Assembler::new(&graph).assemble_contigs();

    // Write output
    fs::create_dir_all("data")?;
    let out_path = Path::new("data").join(&args.output);
    let mut f = fs::File::create(&out_path)?;
    for (i, c) in contigs.iter().enumerate() {
        writeln!(f, ">contig_{}", i + 1)?;
        writeln!(f, "{}", c)?;
    }

    println!("Wrote {} contigs to {:?}", contigs.len(), out_path);
    Ok(())
}
