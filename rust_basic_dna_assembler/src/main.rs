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
    #[clap(short, long, default_value = "0")]
    threshold: usize,
    #[clap(short = 'd', long, default_value = "5")]
    depth: usize,
    #[clap(short = 'b', long, default_value = "20")]
    bubble_depth: usize,
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
    graph.remove_bubbles(args.bubble_depth);

    eprintln!("Loaded {} reads", reads.len());
    let edge_cnt = graph.edge_count();
    let node_cnt = graph.in_degree.len();
    let zero_indeg: Vec<_> = graph
        .in_degree
        .iter()
        .filter_map(|(&u,&d)| if d==0 { Some(u) } else { None })
        .collect();
    eprintln!(
        "Graph: {} edges, {} nodes; {} sources (in_deg=0)",
        edge_cnt, node_cnt, zero_indeg.len()
    );

for &u in zero_indeg.iter().take(5) {
    eprintln!(" source kmer: {}", graph.decode(u));
}

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
