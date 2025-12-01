//! Interfejs wiersza poleceń dla assemblera AI DNA.
//!
//! Ten moduł odpowiada za pełny pipeline składania genomu na podstawie
//! dwustronnych odczytów FASTQ.  
//!
//! Główne etapy działania:
//! - wczytywanie odczytów z plików FASTQ,
//! - budowa grafu de Bruijna na podstawie k-merów,
//! - czyszczenie grafu (filtrowanie, usuwanie zakończeń, usuwanie bąbli),
//! - składanie kontigów z uproszczonego grafu,
//! - zapis wyników do pliku FASTA lub GFA.
//!
//! Obecna wersja zapisuje grafik do GFA,
//! a kod odpowiedzialny za składanie kontigów może być aktywowany w przyszłych krokach.

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

/// Struktura opisująca argumenty wiersza poleceń.
///
/// Opcje umożliwiają konfigurację:
/// - wejściowych plików FASTQ (pierwszy i drugi koniec),
/// - długości k-mera,
/// - progu filtrowania krawędzi,
/// - głębokości usuwania zakończeń,
/// - głębokości usuwania bąbli,
/// - nazwy pliku wynikowego.
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
}

/// Główna funkcja pipeline'u assemblera.
///
/// Odpowiada za:
/// - wczytanie odczytów FASTQ,
/// - zbudowanie grafu k-merowego,
/// - uruchomienie operacji czyszczących,
/// - zapis wyników.
/// 
/// # Zwraca
/// Zwraca `io::Result<()>` — błąd tylko w wypadku problemów z I/O.
fn main() -> io::Result<()> {
    let args = Args::parse();

    // Wczytanie odczytów z dwóch plików FASTQ
    let mut reads = FastqParser::new(&args.reads1).load_reads()?;
    reads.extend(FastqParser::new(&args.reads2).load_reads()?);

    // Utworzenie grafu k-merów
    let mut graph = KmerGraph::new(args.kmer);

    // Budowa grafu de Bruijna
    graph.build(&reads);

    // Etap filtracji i czyszczenia grafu — obecnie wyłączony,
    // ale może zostać aktywowany podczas testów.
    //
    // graph.filter_low_coverage(args.threshold);
    // graph.remove_dead_ends(Some(args.depth));
    // graph.remove_bubbles(args.bubble_depth);

    // Zapis grafu do GFA — aktualnie główny cel etapu Rust
    graph.write_gfa(&args.output).unwrap();

    // Podstawowe statystyki grafu
    eprintln!("Wczytano {} odczytów", reads.len());
    let edge_cnt = graph.edge_count();
    let node_cnt = graph.in_degree.len();

    // Węzły bez poprzedników (źródła)
    let zero_indeg: Vec<_> = graph
        .in_degree
        .iter()
        .filter_map(|(&u, &d)| if d == 0 { Some(u) } else { None })
        .collect();

    eprintln!(
        "Graf: {} krawędzi, {} węzłów; {} źródeł (in_deg=0)",
        edge_cnt, node_cnt, zero_indeg.len()
    );

    // Debug: można wypisać przykładowe k-mery źródłowe
    // for &u in zero_indeg.iter().take(5) {
    //     eprintln!(" źródło: {}", graph.decode(u));
    // }

    // Składanie kontigów — aktualnie wyłączone
    //
    // let contigs = Assembler::new(&graph).assemble_contigs();
    //
    // fs::create_dir_all("data")?;
    // let out_path = Path::new("data").join(&args.output);
    // let mut f = fs::File::create(&out_path)?;
    // for (i, c) in contigs.iter().enumerate() {
    //     writeln!(f, ">contig_{}", i + 1)?;
    //     writeln!(f, "{}", c)?;
    // }
    // println!("Zapisano {} kontigów do {:?}", contigs.len(), out_path);

    Ok(())
}
