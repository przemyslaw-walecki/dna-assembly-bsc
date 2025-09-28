//! # High-k Unitig Scan Labeling (Aho–Corasick over Unitigs)
//!
//! This tool labels BubbleGun bubbles by scanning **high-k** unitig sequences
//! for exact matches to **low-k** candidate path sequences (both strands) using
//! Aho–Corasick.
//!
//! ## Pipeline
//! 1) Read low-k GFA, enumerate orientation-safe candidate sequences per bubble.
//! 2) Read high-k GFA and compact it to unitigs (chains between non-1/1-degree nodes).
//! 3) Build an AC automaton over all candidates (forward+RC), scan all unitigs.
//! 4) Label bubbles with a single winning candidate.
//!
//! ## Winner rule
//! - Default: winner unique and `hits == 1` (`--require-unique-occurrence`).
//! - `--allow-multi-occurrence`: winner unique with `hits ≥ 1`.
//!
//! ## Usage
//! ```text
//! cargo run --release --bin project_highk_unitigs -- \
//!   --low-gfa ecoli_21_graph.gfa \
//!   --bubbles ecoli_unlabeled_bubblegun.json \
//!   --high-gfa ecoli_41_graph.gfa \
//!   --out labels_from_k41.jsonl \
//!   [--allow-multi-occurrence]
//! ```

use aho_corasick::{AhoCorasick, AhoCorasickBuilder};
use serde::Deserialize;
use serde_json::Value;
use std::collections::{HashMap, HashSet, VecDeque};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

// -------------------- shared / util --------------------

/// Parse a decimal string into `u128`, panicking on non-numeric input.
fn parse_u128(s: &str) -> u128 {
    s.parse::<u128>().unwrap_or_else(|_| panic!("Non-numeric id: {s}"))
}

/// Reverse-complement a DNA string (A↔T, C↔G; others→N).
fn revcomp(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.bytes().rev() {
        out.push(match c {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        });
    }
    out
}

/// Append the next *k*-mer `next_seq` to `assembled` if it overlaps by `k-1` bp.
/// Falls back to RC of `next_seq` if needed.
fn append_next_kmer(assembled: &mut String, next_seq: &str, k: usize) -> bool {
    if assembled.as_bytes()[(assembled.len() - (k - 1))..] == next_seq.as_bytes()[..(k - 1)] {
        assembled.push(next_seq.as_bytes()[k - 1] as char);
        true
    } else {
        let rc = revcomp(next_seq);
        if assembled.as_bytes()[(assembled.len() - (k - 1))..] == rc.as_bytes()[..(k - 1)] {
            assembled.push(rc.as_bytes()[k - 1] as char);
            true
        } else {
            false
        }
    }
}

// -------------------- low-k graph + bubblegun --------------------

/// Segment from GFA `S` line.
#[derive(Debug)]
struct Segment {
    id: u128,
    seq: String,
}

/// Parse GFA into segments and both successor/predecessor adjacency maps.
fn parse_gfa(path: &str) -> (HashMap<u128, Segment>, HashMap<u128, Vec<u128>>, HashMap<u128, Vec<u128>>) {
    let f = File::open(path).expect("open GFA");
    let rdr = BufReader::new(f);
    let mut segs: HashMap<u128, Segment> = HashMap::new();
    let mut succ: HashMap<u128, Vec<u128>> = HashMap::new();
    let mut pred: HashMap<u128, Vec<u128>> = HashMap::new();

    for line in rdr.lines().flatten() {
        if line.is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        match cols.get(0).copied() {
            Some("S") => {
                if cols.len() < 3 { continue; }
                let id = parse_u128(cols[1]);
                let seq = cols[2].to_string();
                segs.insert(id, Segment { id, seq });
                succ.entry(id).or_default();
                pred.entry(id).or_default();
            }
            Some("L") => {
                if cols.len() < 6 { continue; }
                let from = parse_u128(cols[1]);
                let to   = parse_u128(cols[3]);
                succ.entry(from).or_default().push(to);
                pred.entry(to).or_default().push(from);
                succ.entry(to).or_default();
                pred.entry(from).or_default();
            }
            _ => {}
        }
    }
    (segs, succ, pred)
}

/// Bubble shapes used for BubbleGun input.
#[derive(Deserialize, Debug)]
struct Bubble { id: usize, ends: Vec<String>, inside: Vec<String> }
#[derive(Deserialize, Debug)]
struct BubbleChain { bubbles: Vec<Bubble> }

/// Load bubbles from a BubbleGun JSON file in flexible shapes.
fn collect_bubbles_from_bubblegun(path: &str) -> Vec<Bubble> {
    let f = File::open(path).expect("open BubbleGun JSON");
    let val: Value = serde_json::from_reader(BufReader::new(f)).expect("parse BubbleGun JSON");

    fn looks_like_bubble(v: &Value) -> bool {
        v.get("id").is_some() && v.get("ends").is_some()
    }

    let mut out = Vec::new();
    match &val {
        Value::Array(a) => {
            if !a.is_empty() && a[0].get("bubbles").is_some() {
                for chain in a {
                    if let Some(bs) = chain.get("bubbles") {
                        if let Ok(v) = serde_json::from_value::<Vec<Bubble>>(bs.clone()) {
                            out.extend(v);
                        }
                    }
                }
            } else {
                for b in a {
                    if looks_like_bubble(b) {
                        if let Ok(bb) = serde_json::from_value::<Bubble>(b.clone()) {
                            out.push(bb);
                        }
                    }
                }
            }
        }
        Value::Object(map) => {
            if let Some(bs) = map.get("bubbles") {
                if let Ok(v) = serde_json::from_value::<Vec<Bubble>>(bs.clone()) { out.extend(v); }
            }
            for v in map.values() {
                if let Some(bs) = v.get("bubbles") {
                    if let Ok(vb) = serde_json::from_value::<Vec<Bubble>>(bs.clone()) {
                        out.extend(vb);
                    }
                }
            }
        }
        _ => {}
    }
    out
}

/// DFS enumerate orientation-safe path sequences inside the bubble’s induced subgraph.
/// Tries `start→end`, else tries reverse direction and RCs the outputs back to start→end.
fn enumerate_paths_sequences(
    start: u128,
    end: u128,
    nodes_set: &HashSet<u128>,
    succ: &HashMap<u128, Vec<u128>>,
    segs: &HashMap<u128, Segment>,
    k: usize,
    max_paths: usize,
    guard_steps: usize,
) -> Vec<String> {
    let mut paths: Vec<String> = Vec::new();

    let mut next_map: HashMap<u128, Vec<u128>> = HashMap::new();
    for &u in nodes_set {
        let mut v = Vec::new();
        if let Some(nxt) = succ.get(&u) {
            for &w in nxt { if nodes_set.contains(&w) { v.push(w); } }
        }
        if !v.is_empty() { next_map.insert(u, v); }
    }

    let mut stack: Vec<(u128, usize, String)> = Vec::new();
    if let Some(vs) = next_map.get(&start) {
        for &v in vs {
            let s0 = segs.get(&start).unwrap().seq.clone();
            let s1 = segs.get(&v).unwrap().seq.clone();
            let mut seq = s0.clone();
            if append_next_kmer(&mut seq, &s1, k) { stack.push((v, 1, seq)); }
        }
    }

    while let Some((at, steps, seq)) = stack.pop() {
        if paths.len() >= max_paths { break; }
        if steps > guard_steps { continue; }
        if at == end { paths.push(seq.clone()); continue; }
        if let Some(vs) = next_map.get(&at) {
            for &v in vs {
                let s1 = &segs.get(&v).unwrap().seq;
                let mut seq2 = seq.clone();
                if append_next_kmer(&mut seq2, s1, k) {
                    stack.push((v, steps + 1, seq2));
                }
            }
        }
    }

    // If none, try reverse direction and RC the results so they’re start→end
    if paths.is_empty() {
        let mut stack: Vec<(u128, usize, String)> = Vec::new();
        if let Some(vs) = next_map.get(&end) {
            for &v in vs {
                let s0 = segs.get(&end).unwrap().seq.clone();
                let s1 = segs.get(&v).unwrap().seq.clone();
                let mut seq = s0.clone();
                if append_next_kmer(&mut seq, &s1, k) { stack.push((v, 1, seq)); }
            }
        }
        while let Some((at, steps, seq)) = stack.pop() {
            if paths.len() >= max_paths { break; }
            if steps > guard_steps { continue; }
            if at == start { paths.push(seq.clone()); continue; }
            if let Some(vs) = next_map.get(&at) {
                for &v in vs {
                    let s1 = &segs.get(&v).unwrap().seq;
                    let mut seq2 = seq.clone();
                    if append_next_kmer(&mut seq2, s1, k) {
                        stack.push((v, steps + 1, seq2));
                    }
                }
            }
        }
        if !paths.is_empty() { paths = paths.into_iter().map(|s| revcomp(&s)).collect(); }
    }

    paths
}

// -------------------- high-k compaction to unitigs --------------------

/// Compact the high-k graph to unitigs (chains between non-1/1-degree nodes)
/// and return the assembled unitig sequences (orientation-safe).
fn compact_highk_unitigs(
    segs: &HashMap<u128, Segment>,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    k: usize,
) -> Vec<String> {
    // compute degrees
    let mut indeg: HashMap<u128, usize> = HashMap::new();
    let mut outdeg: HashMap<u128, usize> = HashMap::new();
    for (&u, vs) in succ {
        outdeg.insert(u, vs.len());
        indeg.entry(u).or_insert(0);
        for &v in vs {
            *indeg.entry(v).or_insert(0) += 1;
            outdeg.entry(v).or_insert(0);
        }
    }

    // supernodes = nodes where deg != 1 (either side)
    let mut supernodes: Vec<u128> = Vec::new();
    for &n in segs.keys() {
        let in1 = *indeg.get(&n).unwrap_or(&0) == 1;
        let out1 = *outdeg.get(&n).unwrap_or(&0) == 1;
        if !(in1 && out1) { supernodes.push(n); }
    }

    let mut unitigs: Vec<String> = Vec::new();

    for &u in &supernodes {
        if let Some(vs0) = succ.get(&u) {
            for &v0 in vs0 {
                // start a chain at (u -> v0)
                let mut chain: Vec<u128> = vec![u, v0];
                let mut cur = v0;
                while *indeg.get(&cur).unwrap_or(&0) == 1 && *outdeg.get(&cur).unwrap_or(&0) == 1 {
                    if let Some(nx) = succ.get(&cur).and_then(|v| v.first()).copied() {
                        chain.push(nx);
                        cur = nx;
                    } else {
                        break;
                    }
                }
                // assemble chain sequence (orientation-safe)
                if chain.len() >= 2 {
                    let mut seq = segs[&chain[0]].seq.clone();
                    let mut ok = true;
                    for node in chain.iter().skip(1) {
                        let s = &segs[node].seq;
                        if !append_next_kmer(&mut seq, s, k) { ok = false; break; }
                    }
                    if !ok {
                        // retry starting with RC of first
                        let mut seq2 = revcomp(&segs[&chain[0]].seq);
                        let mut ok2 = true;
                        for node in chain.iter().skip(1) {
                            let s = &segs[node].seq;
                            if !append_next_kmer(&mut seq2, s, k) { ok2 = false; break; }
                        }
                        if ok2 { unitigs.push(seq2); }
                    } else {
                        unitigs.push(seq);
                    }
                }
            }
        }
    }
    unitigs
}

// -------------------- main --------------------

/// Parse CLI, build low-k candidates, compact high-k to unitigs, AC-scan unitigs,
/// and write JSONL labels (`bubble_id`, `label_path`, `label_reason`).
fn main() {
    // Usage:
    // cargo run --release --bin project_highk_unitigs -- \
    //   --low-gfa ecoli_21_graph.gfa \
    //   --bubbles ecoli_unlabeled_bubblegun.json \
    //   --high-gfa ecoli_41_graph.gfa \
    //   --out labels_from_k41.jsonl \
    //   [--allow-multi-occurrence]
    //
    // Winner rule: exactly one candidate has >0 hits in high-k unitigs.
    // By default also require hits==1; relax with --allow-multi-occurrence.

    let mut low_gfa   = None;
    let mut bubbles   = None;
    let mut high_gfa  = None;
    let mut out_path  = None;
    let mut require_unique_occ = true;

    let mut args = env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--low-gfa"  => low_gfa  = args.next(),
            "--bubbles"  => bubbles  = args.next(),
            "--high-gfa" => high_gfa = args.next(),
            "--out"      => out_path = args.next(),
            "--require-unique-occurrence" => require_unique_occ = true,
            "--allow-multi-occurrence"    => require_unique_occ = false,
            _ => {
                eprintln!("Unknown arg: {a}");
                eprintln!("Expected: --low-gfa <gfa> --bubbles <json> --high-gfa <gfa> --out <jsonl> [--allow-multi-occurrence]");
                std::process::exit(2);
            }
        }
    }
    let low_gfa  = low_gfa.expect("--low-gfa required");
    let bubbles  = bubbles.expect("--bubbles required");
    let high_gfa = high_gfa.expect("--high-gfa required");
    let out_path = out_path.expect("--out required");

    // low-k: for candidate path enumeration
    let (segs_lo, succ_lo, _pred_lo) = parse_gfa(&low_gfa);
    if segs_lo.is_empty() { eprintln!("empty low-k GFA"); std::process::exit(1); }
    let k = segs_lo.values().next().unwrap().seq.len();

    let bub_list = collect_bubbles_from_bubblegun(&bubbles);
    if bub_list.is_empty() { eprintln!("[info] no bubbles found; nothing to do."); std::process::exit(0); }

    // Build low-k candidate patterns
    let max_paths_per_bubble = 512usize;
    let guard_steps_per_path = 10_000usize;

    let mut patterns: Vec<String> = Vec::new();           // forward + RC
    let mut pat_to_key: Vec<(usize, usize)> = Vec::new(); // base -> (bubble_id, path_idx)
    let mut bubble_to_bases: HashMap<usize, Vec<usize>> = HashMap::new();

    eprintln!("[info] enumerating low-k candidates (k={})", k);
    for b in &bub_list {
        if b.ends.len() != 2 { continue; }
        let start = parse_u128(&b.ends[0]);
        let end   = parse_u128(&b.ends[1]);
        if !segs_lo.contains_key(&start) || !segs_lo.contains_key(&end) { continue; }

        let mut nodes: HashSet<u128> = b.inside.iter().map(|s| parse_u128(s)).collect();
        nodes.insert(start); nodes.insert(end);

        let seqs = enumerate_paths_sequences(
            start, end, &nodes, &succ_lo, &segs_lo, k, max_paths_per_bubble, guard_steps_per_path
        );
        if seqs.is_empty() { continue; }

        let mut bases = Vec::new();
        for (pi, s) in seqs.iter().enumerate() {
            if s.is_empty() { continue; }
            let base = pat_to_key.len();
            pat_to_key.push((b.id, pi));
            patterns.push(s.clone());
            patterns.push(revcomp(s));
            bases.push(base);
        }
        if !bases.is_empty() { bubble_to_bases.insert(b.id, bases); }
    }
    if patterns.is_empty() { eprintln!("[info] no candidates assembled; nothing to do."); std::process::exit(0); }

    // high-k: compact to unitigs and assemble unitig sequences
    eprintln!("[info] compacting high-k GFA to unitigs for search");
    let (segs_hi, succ_hi, pred_hi) = parse_gfa(&high_gfa);
    if segs_hi.is_empty() { eprintln!("empty high-k GFA"); std::process::exit(1); }
    let k_hi = segs_hi.values().next().unwrap().seq.len();
    if k_hi < k {
        eprintln!("[warn] high-k ({}) < low-k ({}). You probably swapped inputs.", k_hi, k);
    }
    let unitigs = compact_highk_unitigs(&segs_hi, &succ_hi, &pred_hi, k_hi);
    let total_bp: usize = unitigs.iter().map(|s| s.len()).sum();
    eprintln!("[info] high-k unitigs: {} (total {} bp)", unitigs.len(), total_bp);

    // AC over candidates, scan all unitig sequences
    eprintln!(
        "[info] building Aho–Corasick for {} candidates ({} patterns incl. RC)",
        pat_to_key.len(), patterns.len()
    );
    let ac = AhoCorasickBuilder::new()
        .ascii_case_insensitive(false)
        .build(&patterns)
        .expect("build AC");

    let mut counts: HashMap<usize, usize> = HashMap::new(); // base pattern idx -> hits
    for seq in &unitigs {
        for m in ac.find_iter(seq) {
            let pid = m.pattern().as_usize();
            let base = pid / 2;
            *counts.entry(base).or_insert(0) += 1;
        }
    }

    // pick winners
    eprintln!("[info] deciding winners");
    let mut out = File::create(&out_path).expect("create out");
    let mut labeled = 0usize;

    for (bid, bases) in bubble_to_bases.iter() {
        let mut winner: Option<(usize, usize)> = None; // (base_idx, hits)
        let mut clash = false;

        for &bidx in bases {
            let h = *counts.get(&bidx).unwrap_or(&0);
            if h > 0 {
                if winner.is_none() { winner = Some((bidx, h)); } else { clash = true; break; }
            }
        }
        if let Some((bidx, hits)) = winner {
            if !clash && (!require_unique_occ || hits == 1) {
                let (_bchk, path_idx) = pat_to_key[bidx];
                let reason = if require_unique_occ { "highk_unitig_unique" } else { "highk_unitig_present" };
                let obj = serde_json::json!({"bubble_id": bid, "label_path": path_idx, "label_reason": reason});
                writeln!(out, "{}", obj).unwrap();
                labeled += 1;
            }
        }
    }
    eprintln!("[done] labeled {} bubbles -> {}", labeled, &out_path);
}
