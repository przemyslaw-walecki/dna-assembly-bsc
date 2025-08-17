// main.rs
//
// Hyperbubble labeling pipeline — GNN-ready, ID-free output (parallel + progress).
// - Compact bubbles to superedges (unitigs) with coverage aggregates.
// - Split unitigs at bubble endpoints.
// - Enumerate paths bidirectionally with orientation-safe sequence assembly.
// - Label via adaptive flank extension + exact KMP on both strands.
// - **IDs are removed** from the output. Upstream/downstream and endpoints are sequences.
// - cov_min / cov_mean prefer EC:i per link; fallback to min(KC:u, KC:v).
// - Parallelized over bubbles using Rayon (batched to keep memory bounded).
// - Shows a global progress bar (processed / total / ETA) using `indicatif`.
//
// Usage (positional args):
//   cargo run --release -- <graph.gfa> <bubbles.json/jsonl> <reference.fasta> <out.jsonl>
//
// JSONL per bubble includes:
//   {
//     bubble_id, k,
//     start_seq, end_seq,
//     upstream: [seq,...], downstream: [seq,...],            // <- preset greedy flanks (N=5)
//     upstream_nodes: [{seq,cov},...], downstream_nodes: [{seq,cov},...],
//     nodes: [{seq,cov},...],                  // supernodes in bubble+context
//     edges: [{source_seq,target_seq,len_nodes,len_bp,cov_min,cov_mean}, ...],
//     paths: [[edge_idx,...], ...],
//     label_path: <usize|null>
//   }
//
// Cargo.toml (add):
//   indicatif = "0.17"
//   rayon = "1.10"
//   serde = { version = "1", features = ["derive"] }
//   serde_json = "1"

use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::{Deserializer, Value};
use std::cmp;
use std::collections::{HashMap, HashSet, VecDeque};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

// ---------- Data models ----------

#[derive(Debug)]
struct Segment {
    id: u128,
    seq: String,
    cov: usize, // KC
}

#[derive(Debug, Clone)]
struct Link {
    from: u128,
    to: u128,
    cov: usize, // EC
}

#[derive(Deserialize)]
struct BubbleChain {
    bubbles: Vec<Bubble>,
}

#[derive(Deserialize)]
struct Bubble {
    id: usize,
    ends: Vec<String>,   // [start,end]
    inside: Vec<String>, // node ids
}

#[derive(Clone, Debug)]
struct SuperEdge {
    from: u128,
    to: u128,
    chain: Vec<u128>, // includes endpoints (NOT serialized)
    cov_min: usize,
    cov_mean: f64,
}

// ---------- Output schema (ID-free) ----------

#[derive(Serialize)]
struct BubbleRecord {
    bubble_id: usize,
    k: usize,

    // endpoints as sequences (no IDs)
    start_seq: String,
    end_seq: String,

    // preset flanks (nearest-first, greedy over graph)
    upstream: Vec<String>,
    downstream: Vec<String>,

    // context nodes with coverage (no IDs)
    upstream_nodes: Vec<NodeFeature>,
    downstream_nodes: Vec<NodeFeature>,

    // supernodes inside bubble+context
    nodes: Vec<NodeFeature>,

    // compact superedges (endpoint sequences only)
    edges: Vec<EdgeFeature>,

    // edge-index paths start->end
    paths: Vec<Vec<usize>>,
    label_path: Option<usize>,
}

#[derive(Serialize)]
struct NodeFeature {
    seq: String,
    cov: usize,
}

#[derive(Serialize)]
struct EdgeFeature {
    source_seq: String,
    target_seq: String,
    len_nodes: usize,
    len_bp: usize,
    cov_min: usize,
    cov_mean: f64,
}

// ---------- Utilities ----------

fn parse_u128(s: &str) -> u128 {
    s.parse::<u128>()
        .unwrap_or_else(|_| panic!("Expected numeric node id, got: {s}"))
}

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

// ---------- GFA parsing ----------

fn parse_gfa(
    path: &str,
) -> (HashMap<u128, Segment>, Vec<Link>, HashMap<(u128, u128), usize>) {
    let f = File::open(path).expect("Cannot open GFA");
    let rdr = BufReader::new(f);

    let mut segs: HashMap<u128, Segment> = HashMap::new();
    let mut links: Vec<Link> = Vec::new();
    let mut edge_cov: HashMap<(u128, u128), usize> = HashMap::new();

    for line in rdr.lines().flatten() {
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        match cols.get(0).copied() {
            Some("S") => {
                if cols.len() < 3 {
                    continue;
                }
                let id: u128 = cols[1]
                    .parse()
                    .unwrap_or_else(|_| panic!("Non-numeric segment id: {}", cols[1]));
                let seq = cols[2].to_string();
                let cov = cols.iter()
                .find_map(|f| {
                    let lf = f.to_ascii_lowercase();
                    lf.strip_prefix("kc:i:").and_then(|x| x.parse::<usize>().ok())
                })
                .unwrap_or(0);
                segs.insert(id, Segment { id, seq, cov });
            }
            Some("L") => {
                if cols.len() < 6 {
                    continue;
                }
                let from: u128 = cols[1]
                    .parse()
                    .unwrap_or_else(|_| panic!("Non-numeric link from: {}", cols[1]));
                let to: u128 = cols[3]
                    .parse()
                    .unwrap_or_else(|_| panic!("Non-numeric link to: {}", cols[3]));
                let cov = cols.iter()
                    .find_map(|f| {
                        let lf = f.to_ascii_lowercase();
                        lf.strip_prefix("ec:i:").and_then(|x| x.parse::<usize>().ok())
                    })
                    .unwrap_or(0);
                links.push(Link { from, to, cov });
                edge_cov.insert((from, to), cov);
            }
            _ => {}
        }
    }
    (segs, links, edge_cov)
}

fn build_adj(
    segs: &HashMap<u128, Segment>,
    links: &[Link],
) -> (HashMap<u128, Vec<u128>>, HashMap<u128, Vec<u128>>) {
    let mut succ: HashMap<u128, Vec<u128>> = segs.keys().map(|&k| (k, Vec::new())).collect();
    let mut pred: HashMap<u128, Vec<u128>> = segs.keys().map(|&k| (k, Vec::new())).collect();
    for l in links {
        if let Some(v) = succ.get_mut(&l.from) {
            v.push(l.to);
        }
        if let Some(v) = pred.get_mut(&l.to) {
            v.push(l.from);
        }
    }
    (succ, pred)
}

// ---------- Compaction (coverage aggregation + KC fallback) ----------

fn compact_with_chains(
    nodes: &HashSet<u128>,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    edge_cov: &HashMap<(u128, u128), usize>,
    segments: &HashMap<u128, Segment>,
) -> (Vec<SuperEdge>, Vec<u128>) {
    let mut indeg: HashMap<u128, usize> = HashMap::new();
    let mut outdeg: HashMap<u128, usize> = HashMap::new();
    for &n in nodes {
        indeg.insert(
            n,
            pred.get(&n)
                .map(|v| v.iter().filter(|x| nodes.contains(x)).count())
                .unwrap_or(0),
        );
        outdeg.insert(
            n,
            succ.get(&n)
                .map(|v| v.iter().filter(|x| nodes.contains(x)).count())
                .unwrap_or(0),
        );
    }

    let mut link_cov = |u: u128, v: u128| -> usize {
        let ec = edge_cov.get(&(u, v)).copied().unwrap_or(0);
        if ec > 0 {
            ec
        } else {
            let cu = segments.get(&u).map(|s| s.cov).unwrap_or(0);
            let cv = segments.get(&v).map(|s| s.cov).unwrap_or(0);
            cmp::min(cu, cv)
        }
    };

    let supernodes: Vec<u128> = nodes
        .iter()
        .copied()
        .filter(|n| indeg[n] != 1 || outdeg[n] != 1)
        .collect();

    let mut out = Vec::<SuperEdge>::new();
    for &u in &supernodes {
        if let Some(vs0) = succ.get(&u) {
            for &v0 in vs0 {
                if !nodes.contains(&v0) {
                    continue;
                }
                let mut chain = vec![u, v0];
                let mut cur = v0;

                let mut cov_sum: usize = link_cov(u, v0);
                let mut cov_min: usize = link_cov(u, v0);
                let mut edges_cnt = 1usize;

                while nodes.contains(&cur) && indeg[&cur] == 1 && outdeg[&cur] == 1 {
                    let next = succ.get(&cur).and_then(|v| v.first()).copied();
                    if let Some(nxt) = next {
                        if !nodes.contains(&nxt) {
                            break;
                        }
                        let ec = link_cov(cur, nxt);
                        cov_min = cmp::min(cov_min, ec);
                        cov_sum += ec;
                        edges_cnt += 1;
                        chain.push(nxt);
                        cur = nxt;
                        if indeg[&cur] != 1 || outdeg[&cur] != 1 {
                            break;
                        }
                    } else {
                        break;
                    }
                }
                out.push(SuperEdge {
                    from: u,
                    to: cur,
                    chain,
                    cov_min,
                    cov_mean: if edges_cnt == 0 { 0.0 } else { cov_sum as f64 / edges_cnt as f64 },
                });
            }
        }
    }
    (out, supernodes)
}

// ---------- Splitting at endpoints ----------

fn index_in_chain(chain: &[u128], node: u128) -> Option<usize> {
    chain.iter().position(|&x| x == node)
}

fn split_edges_at(node: u128, edges: &[SuperEdge]) -> Vec<SuperEdge> {
    let mut out = Vec::with_capacity(edges.len() + 2);
    for e in edges {
        if let Some(i) = index_in_chain(&e.chain, node) {
            if i == 0 || i == e.chain.len() - 1 {
                out.push(e.clone());
            } else {
                let left_chain = e.chain[..=i].to_vec();
                let right_chain = e.chain[i..].to_vec();
                out.push(SuperEdge {
                    from: left_chain.first().copied().unwrap(),
                    to: node,
                    chain: left_chain,
                    cov_min: e.cov_min,
                    cov_mean: e.cov_mean,
                });
                out.push(SuperEdge {
                    from: node,
                    to: right_chain.last().copied().unwrap(),
                    chain: right_chain,
                    cov_min: e.cov_min,
                    cov_mean: e.cov_mean,
                });
            }
        } else {
            out.push(e.clone());
        }
    }
    out
}

fn from_index(edges: &[SuperEdge]) -> HashMap<u128, Vec<usize>> {
    let mut m: HashMap<u128, Vec<usize>> = HashMap::new();
    for (i, e) in edges.iter().enumerate() {
        m.entry(e.from).or_default().push(i);
    }
    m
}

// ---------- Orientation-safe assembly ----------

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

fn assemble_chain_sequence(
    chain: &[u128],
    segments: &HashMap<u128, Segment>,
    k: usize,
) -> Option<String> {
    if chain.is_empty() {
        return None;
    }
    // try direct
    let mut seq = segments[&chain[0]].seq.clone();
    let mut ok = true;
    for node in chain.iter().skip(1) {
        let s = &segments[node].seq;
        if !append_next_kmer(&mut seq, s, k) {
            ok = false;
            break;
        }
    }
    if ok {
        return Some(seq);
    }
    // retry with RC of first
    let mut seq2 = revcomp(&segments[&chain[0]].seq);
    for node in chain.iter().skip(1) {
        let s = &segments[node].seq;
        if !append_next_kmer(&mut seq2, s, k) {
            return None;
        }
    }
    Some(seq2)
}

fn append_chain_to_sequence(
    assembled: &mut String,
    chain: &[u128],
    segments: &HashMap<u128, Segment>,
    k: usize,
) -> bool {
    for node in chain.iter().skip(1) {
        let s = &segments[node].seq;
        if !append_next_kmer(assembled, s, k) {
            return false;
        }
    }
    true
}

// ---------- Path enumeration (bidirectional) ----------

fn enumerate_paths_bidir(
    start: u128,
    end: u128,
    edges: &[SuperEdge],
    segments: &HashMap<u128, Segment>,
    k: usize,
    max_paths: usize,
) -> (Vec<Vec<usize>>, Vec<String>, bool) {
    let from_map = from_index(edges);

    fn enumerate_from(
        s: u128,
        t: u128,
        edges: &[SuperEdge],
        from_map: &HashMap<u128, Vec<usize>>,
        segments: &HashMap<u128, Segment>,
        k: usize,
        max_paths: usize,
    ) -> (Vec<Vec<usize>>, Vec<String>) {
        let mut paths: Vec<Vec<usize>> = Vec::new();
        let mut seqs: Vec<String> = Vec::new();
        let mut stack: Vec<(u128, Vec<usize>, String)> = Vec::new();

        if let Some(eidxs) = from_map.get(&s) {
            for &ei in eidxs {
                if let Some(seq) = assemble_chain_sequence(&edges[ei].chain, segments, k) {
                    stack.push((edges[ei].to, vec![ei], seq));
                }
            }
        }

        let guard_max_steps = 2000usize;
        while let Some((at, idxs, seq)) = stack.pop() {
            if paths.len() >= max_paths {
                break;
            }
            if at == t {
                paths.push(idxs.clone());
                seqs.push(seq.clone());
                continue;
            }
            if idxs.len() > guard_max_steps {
                continue;
            }
            if let Some(eidxs) = from_map.get(&at) {
                for &ei in eidxs {
                    let e = &edges[ei];
                    let mut seq2 = seq.clone();
                    if append_chain_to_sequence(&mut seq2, &e.chain, segments, k) {
                        let mut idxs2 = idxs.clone();
                        idxs2.push(ei);
                        stack.push((e.to, idxs2, seq2));
                    }
                }
            }
        }
        (paths, seqs)
    }

    // forward
    let (p_fwd, s_fwd) = enumerate_from(start, end, edges, &from_map, segments, k, max_paths);
    if !p_fwd.is_empty() {
        return (p_fwd, s_fwd, false);
    }
    // reverse
    let (p_rev, s_rev) = enumerate_from(end, start, edges, &from_map, segments, k, max_paths);
    (p_rev, s_rev, true)
}

// ---------- Global degrees + unique flanks (for labeling only) ----------

fn compute_global_deg(
    succ: &HashMap<u128, Vec<u128>>,
    _pred: &HashMap<u128, Vec<u128>>,
) -> (HashMap<u128, usize>, HashMap<u128, usize>) {
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
    (indeg, outdeg)
}

fn gather_unique_chain(
    seed: u128,
    forward: bool, // true successors (right), false predecessors (left)
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    indeg: &HashMap<u128, usize>,
    outdeg: &HashMap<u128, usize>,
    nodes_cap: usize,
) -> Vec<u128> {
    let mut chain: Vec<u128> = vec![seed];
    let mut cur = seed;
    loop {
        if chain.len() >= nodes_cap {
            break;
        }
        let (nexts, cond) = if forward {
            (
                succ.get(&cur).cloned().unwrap_or_default(),
                outdeg.get(&cur).copied().unwrap_or(0) == 1,
            )
        } else {
            (
                pred.get(&cur).cloned().unwrap_or_default(),
                indeg.get(&cur).copied().unwrap_or(0) == 1,
            )
        };
        if !cond || nexts.len() != 1 {
            break;
        }
        let nx = nexts[0];
        chain.push(nx);
        cur = nx;
    }
    chain
}

fn sequence_for_nodes(
    nodes: &[u128],
    segments: &HashMap<u128, Segment>,
    k: usize,
) -> Option<(String, usize)> {
    if nodes.is_empty() {
        return Some(("".to_string(), 0));
    }
    // try forward
    let mut seq = segments[&nodes[0]].seq.clone();
    for node in nodes.iter().skip(1) {
        let s = &segments[node].seq;
        if !append_next_kmer(&mut seq, s, k) {
            // retry with RC of first
            let mut seq2 = revcomp(&segments[&nodes[0]].seq);
            for node2 in nodes.iter().skip(1) {
                let s2 = &segments[node2].seq;
                if !append_next_kmer(&mut seq2, s2, k) {
                    return None;
                }
            }
            let extra2 = if seq2.len() >= k { seq2.len() - k } else { 0 };
            return Some((seq2, extra2));
        }
    }
    let extra = if seq.len() >= k { seq.len() - k } else { 0 };
    Some((seq, extra))
}

fn compute_flanks(
    start: u128,
    end: u128,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    indeg: &HashMap<u128, usize>,
    outdeg: &HashMap<u128, usize>,
    segments: &HashMap<u128, Segment>,
    k: usize,
    max_bp: usize,
) -> (String, String) {
    let left_nodes_rev = gather_unique_chain(
        start,
        /*forward=*/false,
        succ,
        pred,
        indeg,
        outdeg,
        max_bp / cmp::max(k, 1) + 2,
    );
    let left_nodes: Vec<u128> = left_nodes_rev.into_iter().rev().collect();
    let (left_seq_full, _) = sequence_for_nodes(&left_nodes, segments, k).unwrap_or(("".to_string(), 0));
    let right_nodes = gather_unique_chain(
        end,
        /*forward=*/true,
        succ,
        pred,
        indeg,
        outdeg,
        max_bp / cmp::max(k, 1) + 2,
    );
    let (right_seq_full, _) = sequence_for_nodes(&right_nodes, segments, k).unwrap_or(("".to_string(), 0));

    let left_bp = if left_seq_full.len() > k {
        &left_seq_full[..left_seq_full.len() - k]
    } else {
        ""
    };
    let right_bp = if right_seq_full.len() > k {
        &right_seq_full[k..]
    } else {
        ""
    };
    (left_bp.to_string(), right_bp.to_string())
}

// ---------- KMP exact substring search ----------

fn kmp_build(pat: &str) -> Vec<usize> {
    let m = pat.len();
    let pb = pat.as_bytes();
    let mut lps = vec![0usize; m];
    let mut len = 0usize;
    let mut i = 1usize;
    while i < m {
        if pb[i] == pb[len] {
            len += 1;
            lps[i] = len;
            i += 1;
        } else if len != 0 {
            len = lps[len - 1];
        } else {
            lps[i] = 0;
            i += 1;
        }
    }
    lps
}

fn kmp_search_all(pat: &str, text: &str) -> Vec<usize> {
    if pat.is_empty() {
        return vec![];
    }
    let lps = kmp_build(pat);
    let pb = pat.as_bytes();
    let tb = text.as_bytes();
    let mut res = Vec::new();
    let mut i = 0usize;
    let mut j = 0usize;
    while i < tb.len() {
        if tb[i] == pb[j] {
            i += 1;
            j += 1;
            if j == pb.len() {
                res.push(i - j);
                j = lps[j - 1];
            }
        } else if j != 0 {
            j = lps[j - 1];
        } else {
            i += 1;
        }
    }
    res
}

// ---------- Reference text (both strands) ----------

#[derive(Clone)]
struct Chrom {
    name: String,
    seq: String,
    rc_seq: String,
}

struct RefText {
    chroms: Vec<Chrom>,
}

impl RefText {
    fn load_fasta(path: &str) -> Self {
        let f = File::open(path).expect("open FASTA");
        let rdr = BufReader::new(f);
        let mut name = String::new();
        let mut seq = String::new();
        let mut chroms_raw: Vec<(String, String)> = Vec::new();
        for line in rdr.lines().flatten() {
            if line.starts_with('>') {
                if !seq.is_empty() {
                    chroms_raw.push((name.clone(), seq.clone()));
                    seq.clear();
                }
                name = line[1..].trim().to_string();
            } else {
                seq.push_str(line.trim());
            }
        }
        if !seq.is_empty() {
            chroms_raw.push((name, seq));
        }
        let chroms = chroms_raw
            .into_iter()
            .map(|(n, s)| Chrom {
                name: n,
                seq: s.clone(),
                rc_seq: revcomp(&s),
            })
            .collect();
        Self { chroms }
    }
}

// ---------- Labeling with adaptive flanks & KMP ----------

fn label_with_adaptive_flanks(
    candidate_seqs: &[String], // per path (already orientation-safe within bubble)
    start: u128,
    end: u128,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    segments: &HashMap<u128, Segment>,
    k: usize,
    reftext: &RefText,
) -> (Option<usize>, usize) {
    let (indeg, outdeg) = compute_global_deg(succ, pred);
    let max_flank_bp = 4096usize;
    let (left_full, right_full) =
        compute_flanks(start, end, succ, pred, &indeg, &outdeg, segments, k, max_flank_bp);

    let flanks = [0usize, 32, 64, 128, 256, 512, 1024, 2048, 4096];

    for delta in flanks {
        let left_tail = if left_full.len() > delta {
            &left_full[left_full.len() - delta..]
        } else {
            &left_full[..]
        };
        let right_head = if right_full.len() > delta {
            &right_full[..delta]
        } else {
            &right_full[..]
        };

        // total hits per candidate across all chroms (both strands)
        let mut hits_per_path: Vec<usize> = Vec::with_capacity(candidate_seqs.len());
        for seq in candidate_seqs {
            let extended = format!("{left}{mid}{right}", left = left_tail, mid = seq, right = right_head);
            let mut total = 0usize;
            for chrom in reftext.chroms.iter() {
                total += kmp_search_all(&extended, &chrom.seq).len();
                total += kmp_search_all(&extended, &chrom.rc_seq).len();
            }
            hits_per_path.push(total);
        }

        // Unique if exactly one path has exactly one hit, and all others have zero.
        let mut unique_idx: Option<usize> = None;
        let mut bad = false;
        for (i, &h) in hits_per_path.iter().enumerate() {
            match h {
                0 => {}
                1 => {
                    if unique_idx.is_none() { unique_idx = Some(i); }
                    else { bad = true; break; } // two paths with hits -> ambiguous
                }
                _ => { bad = true; break; }   // multi-hits for a path -> ambiguous
            }
        }
        if !bad {
            if let Some(i) = unique_idx {
                return (Some(i), delta * 2);
            }
        }
        // else grow flanks and retry
    }

    (None, max_flank_bp * 2)
}

// ---------- Preset greedy flanks (export only) ----------

fn best_neighbor(
    cur: u128,
    forward: bool,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    segments: &HashMap<u128, Segment>,
    edge_cov: &HashMap<(u128, u128), usize>,
) -> Option<u128> {
    let cand = if forward {
        succ.get(&cur).cloned().unwrap_or_default()
    } else {
        pred.get(&cur).cloned().unwrap_or_default()
    };
    if cand.is_empty() {
        return None;
    }
    // choose by: edge EC (desc), then node KC (desc), then node id (asc) for determinism
    let mut best_id: Option<u128> = None;
    let mut best_key: Option<(usize, usize, u128)> = None; // (ec, kc, id)
    for n in cand {
        let ec = if forward {
            edge_cov.get(&(cur, n)).copied().unwrap_or(0)
        } else {
            edge_cov.get(&(n, cur)).copied().unwrap_or(0)
        };
        let kc = segments.get(&n).map(|s| s.cov).unwrap_or(0);
        let key = (ec, kc, n);
        match best_key {
            None => { best_key = Some(key); best_id = Some(n); }
            Some(bk) => {
                if key.0 > bk.0 || (key.0 == bk.0 && (key.1 > bk.1 || (key.1 == bk.1 && key.2 < bk.2))) {
                    best_key = Some(key);
                    best_id = Some(n);
                }
            }
        }
    }
    best_id
}

fn collect_flanks_greedy(
    seed: u128,
    forward: bool, // upstream: false (use predecessors); downstream: true (use successors)
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    segments: &HashMap<u128, Segment>,
    edge_cov: &HashMap<(u128, u128), usize>,
    n: usize,
) -> Vec<u128> {
    use std::collections::HashSet;
    let mut out: Vec<u128> = Vec::with_capacity(n);
    let mut cur = seed;
    let mut seen: HashSet<u128> = HashSet::new();
    seen.insert(seed);

    for _ in 0..n {
        if let Some(nx) = best_neighbor(cur, forward, succ, pred, segments, edge_cov) {
            if seen.contains(&nx) {
                break; // avoid tiny cycles
            }
            out.push(nx);   // nearest-first
            seen.insert(nx);
            cur = nx;
        } else {
            break;
        }
    }
    out
}

// ---------- Context helpers for JSON ----------

fn collect_nodes_for_bubble_context(
    inside: &HashSet<u128>,
    start: u128,
    end: u128,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    radius: usize,
) -> HashSet<u128> {
    let mut nodes: HashSet<u128> = inside.clone();
    nodes.insert(start);
    nodes.insert(end);

    // upstream expansion
    let mut seen_u: HashSet<u128> = HashSet::new();
    let mut dq = VecDeque::new();
    dq.push_back((start, 0usize));
    seen_u.insert(start);
    while let Some((u, d)) = dq.pop_front() {
        if d == radius {
            continue;
        }
        if let Some(ps) = pred.get(&u) {
            for &p in ps {
                if seen_u.insert(p) {
                    nodes.insert(p);
                    dq.push_back((p, d + 1));
                }
            }
        }
    }

    // downstream expansion
    let mut seen_d: HashSet<u128> = HashSet::new();
    let mut dq2 = VecDeque::new();
    dq2.push_back((end, 0usize));
    seen_d.insert(end);
    while let Some((u, d)) = dq2.pop_front() {
        if d == radius {
            continue;
        }
        if let Some(vs) = succ.get(&u) {
            for &v in vs {
                if seen_d.insert(v) {
                    nodes.insert(v);
                    dq2.push_back((v, d + 1));
                }
            }
        }
    }
    nodes
}

// ---------- Progress helpers ----------

fn count_bubbles_file(path: &str) -> usize {
    let in_file = File::open(path).expect("open bubbles json/jsonl");
    let reader = BufReader::new(in_file);
    let stream = Deserializer::from_reader(reader).into_iter::<Value>();
    let mut total = 0usize;
    for val in stream {
        let Ok(v) = val else { continue }; // skip bad lines
        if v.is_object() {
            let map = v.as_object().unwrap();
            for vv in map.values() {
                if let Ok(bc) = serde_json::from_value::<BubbleChain>(vv.clone()) {
                    total += bc.bubbles.len();
                }
            }
        } else if v.is_array() {
            if let Ok(bcs) = serde_json::from_value::<Vec<BubbleChain>>(v.clone()) {
                for bc in bcs { total += bc.bubbles.len(); }
            } else if let Ok(bs) = serde_json::from_value::<Vec<Bubble>>(v.clone()) {
                total += bs.len();
            }
        } else if let Ok(bc) = serde_json::from_value::<BubbleChain>(v.clone()) {
            total += bc.bubbles.len();
        }
    }
    total
}

// ---------- Per-bubble worker (returns one JSON line) ----------

fn build_record_json(
    bubble: &Bubble,
    k: usize,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    segments: &HashMap<u128, Segment>,
    edge_cov: &HashMap<(u128, u128), usize>,
    reftext: &RefText,
    radius_context: usize,
    flank_kmers: usize,     // <- preset number of kmers to export per side
    max_paths: usize,
) -> Option<String> {
    if bubble.ends.len() != 2 {
        return None;
    }
    let start: u128 = parse_u128(&bubble.ends[0]);
    let end: u128 = parse_u128(&bubble.ends[1]);

    // Induced subgraph nodes (bubble + context)
    let inside: HashSet<u128> = bubble.inside.iter().map(|s| parse_u128(s)).collect();
    let nodes: HashSet<u128> =
        collect_nodes_for_bubble_context(&inside, start, end, succ, pred, radius_context);

    // Compact (with EC fallback to KC)
    let (superedges0, supernodes0) = compact_with_chains(&nodes, succ, pred, edge_cov, segments);

    // Split at endpoints
    let mut superedges = split_edges_at(start, &superedges0);
    superedges = split_edges_at(end, &superedges);

    // Enumerate paths + sequences (bidirectional)
    let (paths_vec, seqs_vec, _reversed) =
        enumerate_paths_bidir(start, end, &superedges, segments, k, max_paths);

    // Label with adaptive flanks + KMP
    let (label_opt, _flank_used_bp) = if paths_vec.is_empty() {
        (None, 0usize)
    } else {
        label_with_adaptive_flanks(&seqs_vec, start, end, succ, pred, segments, k, reftext)
    };

    // ------- Preset flanks to export: greedily take exactly N neighbors if available -------
    let upstream_ids: Vec<u128> = collect_flanks_greedy(
        start, /*forward=*/false, succ, pred, segments, edge_cov, flank_kmers
    );
    let downstream_ids: Vec<u128> = collect_flanks_greedy(
        end, /*forward=*/true, succ, pred, segments, edge_cov, flank_kmers
    );

    // Build JSON record (ID-free)
    let mut rec = BubbleRecord {
        bubble_id: bubble.id, // keep BubbleGun id for bookkeeping
        k,

        start_seq: segments.get(&start).map(|s| s.seq.clone()).unwrap_or_default(),
        end_seq:   segments.get(&end).map(|s| s.seq.clone()).unwrap_or_default(),

        upstream: upstream_ids
            .iter()
            .filter_map(|nid| segments.get(nid).map(|s| s.seq.clone()))
            .collect(),
        downstream: downstream_ids
            .iter()
            .filter_map(|nid| segments.get(nid).map(|s| s.seq.clone()))
            .collect(),

        upstream_nodes: upstream_ids
            .iter()
            .filter_map(|nid| segments.get(nid))
            .map(|s| NodeFeature { seq: s.seq.clone(), cov: s.cov })
            .collect(),
        downstream_nodes: downstream_ids
            .iter()
            .filter_map(|nid| segments.get(nid))
            .map(|s| NodeFeature { seq: s.seq.clone(), cov: s.cov })
            .collect(),

        nodes: Vec::new(),
        edges: Vec::new(),
        paths: paths_vec.clone(),
        label_path: label_opt
    };

    // Supernodes (as sequences)
    for &nid in &supernodes0 {
        if let Some(seg) = segments.get(&nid) {
            rec.nodes.push(NodeFeature {
                seq: seg.seq.clone(),
                cov: seg.cov,
            });
        }
    }

    // Serialize edges with endpoint sequences only
    for e in superedges.iter() {
        let len_nodes = e.chain.len();
        let len_bp = k + len_nodes.saturating_sub(1);

        let source_seq = segments.get(&e.from).map(|s| s.seq.clone()).unwrap_or_default();
        let target_seq = segments.get(&e.to).map(|s| s.seq.clone()).unwrap_or_default();

        rec.edges.push(EdgeFeature {
            source_seq,
            target_seq,
            len_nodes,
            len_bp,
            cov_min: e.cov_min,
            cov_mean: e.cov_mean,
        });
    }

    Some(serde_json::to_string(&rec).unwrap())
}

// ---------- Main (parallel, batched writer, with progress) ----------

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!(
            "Usage: {} <graph.gfa> <bubbles.jsonl> <reference.fasta> <out.jsonl>",
            args.get(0).map(String::as_str).unwrap_or("pipeline")
        );
        std::process::exit(1);
    }
    let gfa_path = &args[1];
    let bubble_path = &args[2];
    let ref_fa = &args[3];
    let out_path = &args[4];

    // Tunables
    let radius_context: usize = 5;   // nodes to include around bubble for compaction
    let flank_kmers: usize = 5;      // <- preset number of upstream/downstream kmers to export
    let max_paths: usize = 256;      // candidate cap per bubble
    const BATCH: usize = 1000;       // bubbles per parallel batch (tweak for memory/CPU)

    // Load graph + reference
    let (segments, links, edge_cov) = parse_gfa(gfa_path);
    if segments.is_empty() {
        eprintln!("Graph has no segments.");
        std::process::exit(1);
    }
    let (succ, pred) = build_adj(&segments, &links);
    let reftext = RefText::load_fasta(ref_fa);

    // k from first segment
    let k = segments.values().next().expect("empty graph").seq.len();
    if k == 0 {
        eprintln!("Detected k=0 from segments; invalid graph.");
        std::process::exit(1);
    }

    // Determine total bubbles for progress (first, quick counting pass)
    let total_bubbles = count_bubbles_file(bubble_path);
    let pb = ProgressBar::new(total_bubbles as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("##-"),
    );

    // IO
    let out_file = File::create(out_path).expect("create output");
    let mut w = BufWriter::new(out_file);
    let in_file = File::open(bubble_path).expect("open bubbles json/jsonl");
    let reader = BufReader::new(in_file);

    // Collect bubbles in streaming fashion, process in parallel batches
    let mut batch: Vec<Bubble> = Vec::with_capacity(BATCH);
    let mut total_written: usize = 0;

    let stream = Deserializer::from_reader(reader).into_iter::<Value>();
    for val in stream {
        let v = match val {
            Ok(v) => v,
            Err(_) => continue,
        };

        // Normalize to list<Bubble>
        let mut bubbles: Vec<Bubble> = Vec::new();
        if v.is_object() {
            let map = v.as_object().unwrap();
            for vv in map.values() {
                let bc: BubbleChain = serde_json::from_value(vv.clone()).expect("cannot parse BubbleChain");
                bubbles.extend(bc.bubbles);
            }
        } else if v.is_array() {
            if let Ok(bcs) = serde_json::from_value::<Vec<BubbleChain>>(v.clone()) {
                for bc in bcs {
                    bubbles.extend(bc.bubbles);
                }
            } else if let Ok(bs) = serde_json::from_value::<Vec<Bubble>>(v.clone()) {
                bubbles.extend(bs);
            } else {
                panic!("Unsupported array JSON structure for bubbles");
            }
        } else if let Ok(bc) = serde_json::from_value::<BubbleChain>(v.clone()) {
            bubbles.extend(bc.bubbles);
        } else {
            // Unknown shape, skip
            continue;
        }

        for b in bubbles {
            batch.push(b);
            if batch.len() == BATCH {
                // Process a batch in parallel; keep per-batch order deterministic
                let pb_for_threads = pb.clone();
                let mut lines: Vec<(usize, String)> = batch
                    .par_iter()
                    .enumerate()
                    .map_init(|| pb_for_threads.clone(), |pbl, (i, bub)| {
                        let out = build_record_json(
                            bub, k, &succ, &pred, &segments, &edge_cov, &reftext,
                            radius_context, flank_kmers, max_paths,
                        );
                        pbl.inc(1); // progress: one bubble done (success or skip)
                        (i, out)
                    })
                    .filter_map(|(i, out)| out.map(|s| (i, s)))
                    .collect();

                lines.sort_by_key(|(i, _)| *i);
                for (_, line) in lines {
                    w.write_all(line.as_bytes()).unwrap();
                    w.write_all(b"\n").unwrap();
                }
                total_written += batch.len();
                batch.clear();
            }
        }
    }

    // flush tail
    if !batch.is_empty() {
        let pb_for_threads = pb.clone();
        let mut lines: Vec<(usize, String)> = batch
            .par_iter()
            .enumerate()
            .map_init(|| pb_for_threads.clone(), |pbl, (i, bub)| {
                let out = build_record_json(
                    bub, k, &succ, &pred, &segments, &edge_cov, &reftext,
                    radius_context, flank_kmers, max_paths,
                );
                pbl.inc(1);
                (i, out)
            })
            .filter_map(|(i, out)| out.map(|s| (i, s)))
            .collect();
        lines.sort_by_key(|(i, _)| *i);
        for (_, line) in lines {
            w.write_all(line.as_bytes()).unwrap();
            w.write_all(b"\n").unwrap();
        }
        total_written += batch.len();
    }

    pb.finish_with_message(format!("Processed {} bubbles", total_bubbles));
    eprintln!("Done: wrote {} hyperbubble records", total_written);
}
