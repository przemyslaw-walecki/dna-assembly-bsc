use aho_corasick::{AhoCorasick, AhoCorasickBuilder};
use serde::{Deserialize};
use serde_json::Value;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::env;

// -------------------- low-k graph --------------------

#[derive(Debug)]
struct Segment {
    id: u128,
    seq: String,
}

fn parse_u128(s: &str) -> u128 {
    s.parse::<u128>().unwrap_or_else(|_| panic!("Non-numeric id: {s}"))
}

fn parse_gfa(path: &str) -> (HashMap<u128, Segment>, HashMap<u128, Vec<u128>>) {
    let f = File::open(path).expect("open GFA");
    let rdr = BufReader::new(f);
    let mut segs: HashMap<u128, Segment> = HashMap::new();
    let mut succ: HashMap<u128, Vec<u128>> = HashMap::new();

    for line in rdr.lines().flatten() {
        if line.is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        match cols.get(0).copied() {
            Some("S") => {
                if cols.len() < 3 { continue; }
                let id: u128 = parse_u128(cols[1]);
                let seq = cols[2].to_string();
                segs.insert(id, Segment { id, seq });
                succ.entry(id).or_default();
            }
            Some("L") => {
                if cols.len() < 6 { continue; }
                let from = parse_u128(cols[1]);
                let to   = parse_u128(cols[3]);
                succ.entry(from).or_default().push(to);
                succ.entry(to).or_default();
            }
            _ => {}
        }
    }
    (segs, succ)
}

// -------------------- BubbleGun parsing --------------------

#[derive(Deserialize, Debug)]
struct Bubble {
    id: usize,
    ends: Vec<String>,     // [start,end]
    inside: Vec<String>,
}

#[derive(Deserialize, Debug)]
struct BubbleChain {
    bubbles: Vec<Bubble>,
}

// Accepts a BubbleGun JSON in any of the common shapes and yields bubbles.
fn collect_bubbles_from_bubblegun(path: &str) -> Vec<Bubble> {
    let f = File::open(path).expect("open BubbleGun JSON");
    let rdr = BufReader::new(f);
    let val: Value = serde_json::from_reader(rdr).expect("parse BubbleGun JSON");

    fn looks_like_bubble(v: &Value) -> bool {
        v.get("id").is_some() && v.get("ends").is_some()
    }

    let mut out: Vec<Bubble> = Vec::new();

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
                // flat list of bubbles
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
            // direct bubbles?
            if let Some(bs) = map.get("bubbles") {
                if let Ok(v) = serde_json::from_value::<Vec<Bubble>>(bs.clone()) {
                    out.extend(v);
                }
            }
            // dict of chains
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

// -------------------- orientation-safe assembly --------------------

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

// Enumerate path sequences within the induced subgraph of {inside ∪ ends}.
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

    // adjacency restricted to nodes_set
    let mut next_map: HashMap<u128, Vec<u128>> = HashMap::new();
    for &u in nodes_set {
        let mut v = Vec::new();
        if let Some(nxt) = succ.get(&u) {
            for &w in nxt {
                if nodes_set.contains(&w) {
                    v.push(w);
                }
            }
        }
        if !v.is_empty() {
            next_map.insert(u, v);
        }
    }

    // seed sequences from start's outgoing edges
    let mut stack: Vec<(u128, usize, String)> = Vec::new(); // (node, depth, seq)
    if let Some(vs) = next_map.get(&start) {
        for &v in vs {
            let s0 = segs.get(&start).unwrap().seq.clone();
            let s1 = segs.get(&v).unwrap().seq.clone();
            let mut seq = s0.clone();
            if append_next_kmer(&mut seq, &s1, k) {
                stack.push((v, 1, seq));
            }
        }
    }

    while let Some((at, steps, seq)) = stack.pop() {
        if paths.len() >= max_paths { break; }
        if steps > guard_steps { continue; }
        if at == end {
            paths.push(seq.clone());
            continue;
        }
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

    // If no path found in this direction, try reverse direction (end→start)
    if paths.is_empty() {
        let mut stack: Vec<(u128, usize, String)> = Vec::new();
        if let Some(vs) = next_map.get(&end) {
            for &v in vs {
                let s0 = segs.get(&end).unwrap().seq.clone();
                let s1 = segs.get(&v).unwrap().seq.clone();
                let mut seq = s0.clone();
                if append_next_kmer(&mut seq, &s1, k) {
                    stack.push((v, 1, seq));
                }
            }
        }
        while let Some((at, steps, seq)) = stack.pop() {
            if paths.len() >= max_paths { break; }
            if steps > guard_steps { continue; }
            if at == start {
                paths.push(seq.clone());
                continue;
            }
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
        // reverse-complement assembled sequences so they correspond to start→end
        if !paths.is_empty() {
            paths = paths.into_iter().map(|s| revcomp(&s)).collect();
        }
    }

    paths
}

// -------------------- high-k scan --------------------

fn stream_gfa_and_count_matches(
    gfa_path: &str,
    ac: &AhoCorasick,
    counts: &mut HashMap<usize, usize>, // by base pattern index
) {
    let f = File::open(gfa_path).expect("open high-k GFA");
    let rdr = BufReader::new(f);
    for line in rdr.lines().flatten() {
        if !line.starts_with('S') { continue; }
        let mut it = line.split('\t');
        if it.next() != Some("S") { continue; }
        let _sid = it.next();
        let seq = it.next();
        if let Some(s) = seq {
            for m in ac.find_iter(s) {
                let pid = m.pattern().as_usize(); // even=forward, odd=RC
                let base = pid / 2;
                *counts.entry(base).or_insert(0) += 1;
            }
        }
    }
}

// -------------------- main --------------------

fn main() {
    // Usage:
    // cargo run --release --bin project_highk_from_bubblegun -- \
    //   --low-gfa k21.gfa \
    //   --bubbles k21_unlabeled_bubblegun.json \
    //   --high-gfa k41.gfa \
    //   --out labels_from_k41.jsonl \
    //   [--allow-multi-occurrence]
    //
    // Labels a bubble when exactly one candidate path has hits in high-k GFA.
    // By default, also require winner hits==1; relax with --allow-multi-occurrence.

    let mut low_gfa   = None;
    let mut bubbles   = None;
    let mut high_gfa  = None;
    let mut out_path  = None;
    let mut require_unique_occ = true;

    let mut args = env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--low-gfa"   => low_gfa  = args.next(),
            "--bubbles"   => bubbles  = args.next(),
            "--high-gfa"  => high_gfa = args.next(),
            "--out"       => out_path = args.next(),
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

    eprintln!("[info] reading low-k graph: {}", &low_gfa);
    let (segs, succ) = parse_gfa(&low_gfa);
    if segs.is_empty() { eprintln!("empty low-k GFA"); std::process::exit(1); }
    let k = segs.values().next().unwrap().seq.len();
    if k < 2 { eprintln!("invalid k from GFA"); std::process::exit(1); }

    eprintln!("[info] loading unlabeled BubbleGun bubbles: {}", &bubbles);
    let bub_list = collect_bubbles_from_bubblegun(&bubbles);
    if bub_list.is_empty() {
        eprintln!("[info] no bubbles found in BubbleGun file; nothing to do.");
        std::process::exit(0);
    }

    // Build patterns from all candidate paths across all bubbles
    let max_paths_per_bubble = 256usize;
    let guard_steps_per_path = 4_000usize; // cap DFS depth by node steps

    let mut patterns: Vec<String> = Vec::new();           // forward and RC interleaved
    let mut pat_to_key: Vec<(usize, usize)> = Vec::new(); // base pattern idx -> (bubble_id, path_idx)
    let mut bubble_to_bases: HashMap<usize, Vec<usize>> = HashMap::new();

    eprintln!("[info] enumerating candidate paths & assembling sequences (k={})", k);
    for b in bub_list.iter() {
        if b.ends.len() != 2 { continue; }
        let start = parse_u128(&b.ends[0]);
        let end   = parse_u128(&b.ends[1]);

        // Induced nodes: inside ∪ {start,end} that exist in the low-k graph.
        let mut nodes_set: HashSet<u128> =
            b.inside.iter().map(|s| parse_u128(s)).filter(|n| segs.contains_key(n)).collect();
        if !segs.contains_key(&start) || !segs.contains_key(&end) {
            continue;
        }
        nodes_set.insert(start);
        nodes_set.insert(end);

        let seqs = enumerate_paths_sequences(
            start, end, &nodes_set, &succ, &segs, k,
            max_paths_per_bubble, guard_steps_per_path
        );
        if seqs.is_empty() { continue; }

        let mut bases: Vec<usize> = Vec::new();
        for (pi, s) in seqs.iter().enumerate() {
            if s.is_empty() { continue; }
            // base index for forward; RC is base+1 in patterns vector
            let base = pat_to_key.len();
            pat_to_key.push((b.id, pi));
            patterns.push(s.clone());
            patterns.push(revcomp(s));
            bases.push(base);
        }
        if !bases.is_empty() {
            bubble_to_bases.insert(b.id, bases);
        }
    }

    if patterns.is_empty() {
        eprintln!("[info] no candidate sequences assembled; nothing to do.");
        std::process::exit(0);
    }

    eprintln!(
        "[info] building Aho-Corasick over {} candidates ({} patterns incl. RC)",
        pat_to_key.len(),
        patterns.len()
    );
    let ac = AhoCorasickBuilder::new()
        .ascii_case_insensitive(false)
        .build(&patterns)
        .expect("build Aho-Corasick");

    eprintln!("[info] streaming high-k GFA: {}", &high_gfa);
    let mut counts: HashMap<usize, usize> = HashMap::new();
    stream_gfa_and_count_matches(&high_gfa, &ac, &mut counts);

    eprintln!("[info] deciding winners");
    let mut out = File::create(&out_path).expect("create out");
    let mut labeled = 0usize;

    for (bid, bases) in bubble_to_bases.iter() {
        let mut winner: Option<(usize, usize)> = None; // (base_idx, hits)
        let mut clash = false;

        for &bidx in bases {
            let h = *counts.get(&bidx).unwrap_or(&0);
            if h > 0 {
                if winner.is_none() {
                    winner = Some((bidx, h));
                } else {
                    clash = true;
                    break;
                }
            }
        }
        if let Some((bidx, hits)) = winner {
            if !clash && (!require_unique_occ || hits == 1) {
                let (_bubble_check, path_idx) = pat_to_key[bidx];
                let reason = if require_unique_occ { "highk_unitig_unique" } else { "highk_unitig_present" };
                let obj = serde_json::json!({
                    "bubble_id": bid,
                    "label_path": path_idx,
                    "label_reason": reason
                });
                writeln!(out, "{}", obj).unwrap();
                labeled += 1;
            }
        }
    }

    eprintln!("[done] labeled {} bubbles -> {}", labeled, &out_path);
}
