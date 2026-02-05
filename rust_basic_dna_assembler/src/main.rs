//! Interfejs wiersza poleceń dla assemblera AI DNA.
//!
//! This binary now supports an external bubble-resolution pipeline:
//!   1) export GFA
//!   2) run BubbleGun CLI on GFA to produce bubbles JSONL
//!   3) run your Python GNN on that JSONL to produce decisions JSONL
//!   4) import decisions and resolve bubbles in-graph
//!   5) assemble unitigs and write FASTA
//!
//! The legacy bubble remover is left as dead code in `kmer_graph.rs`.

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
use std::process::Command;

/// Command-line arguments for the assembler.
#[derive(Parser)]
struct Args {
    /// Ścieżka do pierwszego pliku FASTQ (`-1`)
    #[clap(short = '1')]
    reads1: String,

    /// Ścieżka do drugiego pliku FASTQ (`-2`)
    #[clap(short = '2')]
    reads2: String,

    /// Długość k-mera wykorzystywana w grafie de Bruijna (`-k`)
    #[clap(short, long, default_value = "55")]
    kmer: usize,

    /// Minimalne pokrycie krawędzi, które pozostaje w grafie (`-t`)
    #[clap(short, long, default_value = "0")]
    threshold: usize,

    /// Maksymalna głębokość usuwania zakończeń (`-d`)
    #[clap(short = 'd', long, default_value = "5")]
    depth: usize,

    /// Maksymalna głębokość usuwania bąbli (`-b`)
    #[clap(short = 'b', long, default_value = "20")]
    bubble_depth: usize,

    /// Nazwa pliku wynikowego (`-o`)
    #[clap(short = 'o', long, default_value = "assembled.fa")]
    output: String,

    // ------------------- External pipeline knobs -------------------

    /// Path to write the working GFA snapshot (e.g., `work/graph.gfa`)
    #[clap(long, default_value = "work/graph.gfa")]
    gfa_out: String,

    /// Optional: path to BubbleGun executable (if provided, we attempt to run it)
    #[clap(long)]
    bubblegun_bin: Option<String>,

    /// BubbleGun output (bubbles JSONL). If not provided, we infer `work/bubbles.jsonl`
    #[clap(long)]
    bubbles_jsonl: Option<String>,

    /// Optional: path to Python interpreter (default: `python`)
    #[clap(long, default_value = "python")]
    python_bin: String,

    /// Optional: path to your GNN runner script (if provided, we attempt to run it)
    #[clap(long)]
    gnn_script: Option<String>,

    /// GNN output decisions JSONL. If not provided, we infer `work/decisions.jsonl`
    #[clap(long)]
    decisions_jsonl: Option<String>,

    /// Pass-through extra args to BubbleGun (quoted as one string).
    #[clap(long)]
    bubblegun_args: Option<String>,

    /// Pass-through extra args to the GNN script (quoted as one string).
    #[clap(long)]
    gnn_args: Option<String>,

    /// If set, skip running BubbleGun (assume `--bubbles-jsonl` already exists)
    #[clap(long)]
    skip_bubblegun: bool,

    /// If set, skip running GNN (assume `--decisions-jsonl` already exists)
    #[clap(long)]
    skip_gnn: bool,

    /// If set, perform dead-end trimming before export
    #[clap(long)]
    trim_tips: bool,

    /// Bubble resolution mode:
    ///   - "gnn": run external Python GNN and apply decisions JSONL
    ///   - "coverage": Rust heuristic (keep highest-coverage edges in each bubble)
    ///   - "none": do not resolve bubbles
    #[clap(long, default_value = "gnn")]
    bubble_resolver: String,
}

fn run_bubblegun(
    bubblegun_bin: &str,
    gfa_in: &str,
    bubbles_out: &str,
    extra_args: Option<&str>,
) -> io::Result<()> {
    let mut cmd = Command::new(bubblegun_bin);

    cmd.arg("-g").arg(gfa_in);
    cmd.arg("bchains");
    cmd.arg("--bubble_json").arg(bubbles_out);

    // Optional pass-through extra args
    if let Some(extra) = extra_args {
        for tok in extra.split_whitespace() {
            cmd.arg(tok);
        }
    }

    eprintln!("[pipeline] Running BubbleGun: {:?}", cmd);
    let status = cmd.status()?;
    if !status.success() {
        Err(io::Error::new(
            io::ErrorKind::Other,
            format!("BubbleGun exited with status {}", status),
        ))
    } else {
        Ok(())
    }
}

fn run_gnn_script(
    python_bin: &str,
    gnn_script: &str,
    bubbles_in: &str,
    decisions_out: &str,
    extra_args: Option<&str>,
) -> io::Result<()> {
    let mut cmd = Command::new(python_bin);
    cmd.arg(gnn_script);
    // Adjust flags to what your runner expects:
    // e.g., python gnn_runner.py --in bubbles.jsonl --out decisions.jsonl
    cmd.arg("--in").arg(bubbles_in);
    cmd.arg("--out").arg(decisions_out);
    if let Some(extra) = extra_args {
        for tok in extra.split_whitespace() {
            cmd.arg(tok);
        }
    }
    let status = cmd.status()?;
    if !status.success() {
        Err(io::Error::new(
            io::ErrorKind::Other,
            format!("GNN script exited with status {}", status),
        ))
    } else {
        Ok(())
    }
}

/// Główna funkcja pipeline'u assemblera.
///
/// Pipeline:
///   - load FASTQs
///   - build graph
///   - optional filtering / tip trimming
///   - write GFA
///   - run BubbleGun (unless skipped)
///   - run GNN (unless skipped)
///   - apply decisions to resolve bubbles
///   - assemble unitigs and write FASTA
fn main() -> io::Result<()> {
    let args = Args::parse();

    // Simple safety: coverage-based resolver requires BubbleGun output.
    if args.bubble_resolver == "coverage" && args.skip_bubblegun {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--bubble-resolver=coverage requires BubbleGun; do not pass --skip-bubblegun",
        ));
    }

    // Load reads
    let mut reads = FastqParser::new(&args.reads1).load_reads()?;
    reads.extend(FastqParser::new(&args.reads2).load_reads()?);

    // Build graph
    let mut graph = KmerGraph::new(args.kmer);
    graph.build(&reads);

    if args.threshold > 0 {
        graph.filter_low_coverage(args.threshold);
    }
    if args.trim_tips {
        graph.remove_dead_ends(Some(args.depth));
    }

    // Ensure work dir exists
    if let Some(dir) = Path::new(&args.gfa_out).parent() {
        if !dir.as_os_str().is_empty() {
            fs::create_dir_all(dir)?;
        }
    }

    // Write GFA for downstream BubbleGun
    graph.write_gfa(&args.gfa_out)?;

    // BubbleGun output path
    let bubbles_jsonl = args
        .bubbles_jsonl
        .unwrap_or_else(|| "work/bubbles.jsonl".to_string());
    if let Some(dir) = Path::new(&bubbles_jsonl).parent() {
        if !dir.as_os_str().is_empty() {
            fs::create_dir_all(dir)?;
        }
    }

    // Run BubbleGun unless skipped
    if !args.skip_bubblegun {
        let Some(bin) = args.bubblegun_bin.as_deref() else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Missing --bubblegun-bin (or pass --skip-bubblegun with an existing --bubbles-jsonl).",
            ));
        };
        eprintln!("[pipeline] running BubbleGun…");
        run_bubblegun(bin, &args.gfa_out, &bubbles_jsonl, args.bubblegun_args.as_deref())?;
        eprintln!("[pipeline] BubbleGun done -> {}", bubbles_jsonl);
    } else {
        eprintln!("[pipeline] skipping BubbleGun; using {}", bubbles_jsonl);
    }

    // GNN decisions path
    // (only used when --bubble-resolver=gnn; other modes skip GNN and use heuristics or no resolution)
    let decisions_jsonl = args
        .decisions_jsonl
        .unwrap_or_else(|| "work/decisions.jsonl".to_string());
    if let Some(dir) = Path::new(&decisions_jsonl).parent() {
        if !dir.as_os_str().is_empty() {
            fs::create_dir_all(dir)?;
        }
    }

    // Choose bubble resolution strategy based on CLI flag.
    match args.bubble_resolver.as_str() {
        "gnn" => {
            // Run GNN unless skipped
            if !args.skip_gnn {
                let Some(script) = args.gnn_script.as_deref() else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Missing --gnn-script (or pass --skip-gnn with an existing --decisions-jsonl).",
                    ));
                };
                eprintln!("[pipeline] running GNN…");
                run_gnn_script(
                    &args.python_bin,
                    script,
                    &bubbles_jsonl,
                    &decisions_jsonl,
                    args.gnn_args.as_deref(),
                )?;
                eprintln!("[pipeline] GNN done -> {}", decisions_jsonl);
            } else {
                eprintln!("[pipeline] skipping GNN; using {}", decisions_jsonl);
            }

            // Apply GNN decisions to the in-memory graph
            eprintln!("[pipeline] applying decisions…");
            let removed = graph.resolve_bubbles_from_jsonl(&decisions_jsonl)?;
            eprintln!("[pipeline] removed {} edges via GNN decisions", removed);
        }

        "coverage" => {
            // Pure Rust heuristic: for each BubbleGun bubble, keep the highest-EC edge per source.
            eprintln!("[pipeline] resolving bubbles via coverage heuristic…");
            let removed = graph.resolve_bubbles_by_coverage_from_bubblegun(&bubbles_jsonl)?;
            eprintln!(
                "[pipeline] removed {} edges via coverage-based heuristic",
                removed
            );
        }

        "none" => {
            // No bubble resolution – you just use the cleaned graph as-is.
            eprintln!("[pipeline] bubble resolution disabled (--bubble-resolver=none)");
        }

        other => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Unknown --bubble-resolver mode: {other} (expected gnn|coverage|none)"),
            ));
        }
    }

    // Assemble unitigs and write FASTA
    let contigs = Assembler::new(&graph).assemble_unitigs();

    if let Some(dir) = Path::new(&args.output).parent() {
        if !dir.as_os_str().is_empty() {
            fs::create_dir_all(dir)?;
        }
    }
    let mut f = fs::File::create(&args.output)?;
    for (i, c) in contigs.iter().enumerate() {
        writeln!(f, ">contig_{}", i + 1)?;
        writeln!(f, "{}", c)?;
    }
    eprintln!("[pipeline] wrote {} contigs to {}", contigs.len(), args.output);

    Ok(())
}
