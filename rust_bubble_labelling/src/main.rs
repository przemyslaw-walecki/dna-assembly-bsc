//! # Hyperbubble Labeling – program główny
//!
//! Ten plik:
//! - parsuje argumenty wiersza poleceń,
//! - wczytuje graf GFA, plik bąbli (JSON/JSONL) oraz sekwencję referencyjną,
//! - konfiguruję parametry pipeline'u,
//! - przetwarza bąble w partiach (batch) równolegle przy użyciu Rayona,
//! - zapisuje wynik w formacie JSONL (1 rekord na bąbel).

mod gfa;
mod reference;
mod bubble;

use crate::gfa::{parse_gfa, build_adj};
use crate::reference::RefText;
use crate::bubble::{Bubble, BubbleChain, build_record_json, count_bubbles_file};

use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde_json::{Deserializer, Value};

use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::time::{Duration, Instant};

/// Limit czasu na przetworzenie pojedynczego bąbla (w sekundach).
const PER_BUBBLE_TIMEOUT_SECS: u64 = 120;

/// Główna funkcja programu.
///
/// Wejście:
/// - `graph.gfa` – skompaktowany graf De Bruijna,
/// - `bubbles.json|jsonl` – opis bąbli,
/// - `reference.fasta` – sekwencja referencyjna,
/// - `out.jsonl` – plik wyjściowy.
///
/// Każdy bąbel jest zamieniany na rekord JSONL zawierający:
/// - sekwencje węzłów i krawędzi,
/// - flanki upstream/downstream,
/// - kandydackie ścieżki przez bąbel,
/// - (opcjonalnie) jedną z nich oznaczoną jako ścieżka referencyjna.
fn main() {
    // ---------------------- Parsowanie argumentów ----------------------
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!(
            "Użycie: {} <graph.gfa> <bubbles.json|jsonl> <reference.fasta> <out.jsonl>",
            args.get(0).map(String::as_str).unwrap_or("hyperbubble-label")
        );
        std::process::exit(1);
    }

    let gfa_path = &args[1];
    let bubble_path = &args[2];
    let ref_fa = &args[3];
    let out_path = &args[4];

    // ---------------------- Parametry pipeline'u ----------------------
    //
    // Można je dostosować w zależności od wielkości danych / dostępnej pamięci.
    let radius_context: usize = 12; // promień (w liczbie węzłów) użyty do zbudowania kontekstu bąbla
    let flank_kmers: usize = 5;     // liczba k-merów w górę i w dół, eksportowana jako flanki
    let max_paths: usize = 2048;    // maksymalna liczba enumerowanych ścieżek na bąbel
    const BATCH: usize = 100;       // liczba bąbli przetwarzanych równolegle w jednym batchu

    // ---------------------- Ładowanie grafu GFA ----------------------
    let (segments, links, edge_cov) = parse_gfa(gfa_path);
    if segments.is_empty() {
        eprintln!("Graf nie zawiera segmentów (plik GFA jest pusty lub niepoprawny).");
        std::process::exit(1);
    }

    // Budowa list następników/poprzedników (sukcesorów i poprzedników).
    let (succ, pred) = build_adj(&segments, &links);

    // ---------------------- Ładowanie referencji ----------------------
    let reftext = RefText::load_fasta(ref_fa);

    // Zakładamy, że długość sekwencji w pierwszym segmencie odpowiada k.
    let k = segments
        .values()
        .next()
        .expect("Pusty graf – brak segmentów")
        .seq
        .len();
    if k == 0 {
        eprintln!("Wykryto k=0 na podstawie segmentów – niepoprawny graf.");
        std::process::exit(1);
    }

    // ---------------------- Pasek postępu ----------------------
    let total_bubbles = count_bubbles_file(bubble_path);
    let pb = ProgressBar::new(total_bubbles as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("##-"),
    );

    // ---------------------- Wyjściowy plik JSONL ----------------------
    let out_file =
        File::create(out_path).expect("Nie można utworzyć pliku wyjściowego JSONL");
    let mut w = BufWriter::new(out_file);

    // ---------------------- Wejściowy plik z bąblami ----------------------
    let in_file = File::open(bubble_path).expect("Nie można otworzyć pliku z bąblami");
    let reader = BufReader::new(in_file);

    // Będziemy czytać strumień JSON (JSON/JSONL) i normalizować go do listy `Bubble`.
    let mut batch: Vec<Bubble> = Vec::with_capacity(BATCH);
    let stream = Deserializer::from_reader(reader).into_iter::<Value>();

    for val in stream {
        let Ok(v) = val else { continue };

        // Normalizacja do Vec<Bubble> – obsługa kilku wariantów formatu:
        // - pojedynczy obiekt zawierający BubbleChain,
        // - tablica BubbleChain,
        // - tablica Bubble.
        let mut bubbles: Vec<Bubble> = Vec::new();

        if v.is_object() {
            let map = v.as_object().unwrap();
            for vv in map.values() {
                if let Ok(bc) = serde_json::from_value::<BubbleChain>(vv.clone()) {
                    bubbles.extend(bc.bubbles);
                }
            }
        } else if v.is_array() {
            if let Ok(bcs) = serde_json::from_value::<Vec<BubbleChain>>(v.clone()) {
                for bc in bcs {
                    bubbles.extend(bc.bubbles);
                }
            } else if let Ok(bs) = serde_json::from_value::<Vec<Bubble>>(v.clone()) {
                bubbles.extend(bs);
            } else {
                // Inny format tablicy – pomijamy.
                continue;
            }
        } else if let Ok(bc) = serde_json::from_value::<BubbleChain>(v.clone()) {
            bubbles.extend(bc.bubbles);
        } else {
            // Nieobsługiwany kształt JSON – pomijamy.
            continue;
        }

        // Dodajemy bąble do batcha i przetwarzamy, gdy batch osiągnie rozmiar BATCH.
        for b in bubbles {
            batch.push(b);
            if batch.len() == BATCH {
                let curr_batch: Vec<Bubble> = std::mem::take(&mut batch);
                let pb_for_threads = pb.clone();

                // Przetwarzanie w Rayonie – każdy bąbel niezależnie, do deadline’u.
                let mut lines: Vec<(usize, String)> = curr_batch
                    .par_iter()
                    .enumerate()
                    .map_init(
                        || pb_for_threads.clone(),
                        |pbl, (i, bub)| {
                            let deadline = Instant::now()
                                + Duration::from_secs(PER_BUBBLE_TIMEOUT_SECS);
                            let out = build_record_json(
                                bub,
                                k,
                                &succ,
                                &pred,
                                &segments,
                                &edge_cov,
                                &reftext,
                                radius_context,
                                flank_kmers,
                                max_paths,
                                deadline,
                            );
                            pbl.inc(1);
                            (i, out)
                        },
                    )
                    .filter_map(|(i, out)| out.map(|s| (i, s)))
                    .collect();

                // Zachowujemy kolejność z batcha (sort po indeksie).
                lines.sort_by_key(|(i, _)| *i);
                for (_, line) in lines {
                    w.write_all(line.as_bytes()).unwrap();
                    w.write_all(b"\n").unwrap();
                }
            }
        }
    }

    // ---------------------- Przetworzenie „ogona” batcha ----------------------
    if !batch.is_empty() {
        let curr_batch: Vec<Bubble> = std::mem::take(&mut batch);
        let pb_for_threads = pb.clone();

        let mut lines: Vec<(usize, String)> = curr_batch
            .par_iter()
            .enumerate()
            .map_init(
                || pb_for_threads.clone(),
                |pbl, (i, bub)| {
                    let deadline =
                        Instant::now() + Duration::from_secs(PER_BUBBLE_TIMEOUT_SECS);
                    let out = build_record_json(
                        bub,
                        k,
                        &succ,
                        &pred,
                        &segments,
                        &edge_cov,
                        &reftext,
                        radius_context,
                        flank_kmers,
                        max_paths,
                        deadline,
                    );
                    pbl.inc(1);
                    (i, out)
                },
            )
            .filter_map(|(i, out)| out.map(|s| (i, s)))
            .collect();

        lines.sort_by_key(|(i, _)| *i);
        for (_, line) in lines {
            w.write_all(line.as_bytes()).unwrap();
            w.write_all(b"\n").unwrap();
        }
    }

    pb.finish_with_message(format!("Przetworzono {} bąbli", total_bubbles));
}
