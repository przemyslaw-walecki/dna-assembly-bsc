//! # Hyperbubble Labeling Pipeline (De Bruijn, GNN-ready, ID-free)
//!
//! This executable takes a compacted De Bruijn graph (GFA), a set of hyperbubbles,
//! and a reference genome (FASTA), then emits **one JSONL record per bubble**
//! with **sequence-only features** suitable for GNN training. Node/edge IDs are
//! removed — sequences are used instead.
//!
//! ## Highlights
//! - Compacts simple chains into **superedges (unitigs)** and aggregates coverage.
//! - **Splits** unitigs at bubble endpoints to preserve boundary semantics.
//! - Enumerates **all paths** between endpoints (bounded), assembling sequences
//!   in an **orientation-safe** manner (both strands).
//! - **Primary labeler:** reference-anchored walk between endpoint *k*-mers (both strands).
//!   - Fallback to **Minimal Unique Extension (MUE)** on the reference when endpoints are repetitive.
//! - **Fallback labeler (lenient, multi-locus, 100% safe):**
//!   - Adaptive flank extension + **Aho–Corasick** multi-pattern scan over forward reference.
//!   - Collapse candidates by reverse-complement class. **Accept** when **exactly one class** has ≥1 hit.
//!   - If the winning class has multiple candidate indices with identical sequence, set
//!     `label_reason="adaptive_kmp_sequence_tie"` and **omit** `label_path` (safety gate).
//! - **Timeout safety:** each bubble has a hard **120s** *cooperative* deadline. If exceeded at any point,
//!   emit a minimal record with `label_reason="timeout"` (no label path).
//! - **ID-free output:** upstream/downstream and endpoints are sequences.
//! - Coverage: `cov_min` / `cov_mean` prefer **EC:i** per link; fallback to `min(KC:u, KC:v)`.
//! - **Parallelized** with Rayon (batched to bound memory); global progress via `indicatif`.
//!
//! ## Usage
//! ```text
//! cargo run --release -- <graph.gfa> <bubbles.json|jsonl> <reference.fasta> <out.jsonl>
//! ```
//!
//! ## Cargo dependency (Cargo.toml)
//! ```toml
//! aho-corasick = "1.1"
//! rayon = "1.10"
//! indicatif = "0.17"
//! serde = { version = "1", features = ["derive"] }
//! serde_json = "1.0"
//! ```

use aho_corasick::{AhoCorasickBuilder, MatchKind};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::{Deserializer, Value};
use std::cmp;
use std::collections::{HashMap, HashSet, VecDeque};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time::{Duration, Instant};

const PER_BUBBLE_TIMEOUT_SECS: u64 = 120;

/// Segment (GFA `S` record).
#[derive(Debug)]
struct Segment {
    id: u128,
    seq: String,
    cov: usize, // KC
}

/// Link (GFA `L` record).
#[derive(Debug, Clone)]
struct Link {
    from: u128,
    to: u128,
    cov: usize, // EC
}

/// Container for a chain of bubbles (input JSON helper).
#[derive(Deserialize)]
struct BubbleChain {
    bubbles: Vec<Bubble>,
}

/// A single hyperbubble definition.
#[derive(Deserialize)]
struct Bubble {
    id: usize,
    ends: Vec<String>,   // [start,end]
    inside: Vec<String>, // node ids
}

/// A compacted path (“superedge”/unitig) across nodes with coverage aggregates.
#[derive(Clone, Debug)]
struct SuperEdge {
    from: u128,
    to: u128,
    chain: Vec<u128>, // includes endpoints (NOT serialized)
    cov_min: usize,
    cov_mean: f64,
}

/// JSONL record: ID-free bubble features suitable for GNN ingestion.
#[derive(Serialize)]
struct BubbleRecord {
    bubble_id: usize,
    k: usize,
    start_seq: String,
    end_seq: String,
    upstream: Vec<String>,
    downstream: Vec<String>,
    upstream_nodes: Vec<NodeFeature>,
    downstream_nodes: Vec<NodeFeature>,
    nodes: Vec<NodeFeature>,
    edges: Vec<EdgeFeature>,
    // NEW: strategy A (meta) inside picks
    inside_nodes: Vec<InsideNodeFeature>,
    paths: Vec<Vec<usize>>,
    label_path: Option<usize>,
    label_reason: Option<String>,
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

// -------- Strategy A: topology-aware inside node selection (no labels/ref) --------

#[derive(Serialize, Clone)]
struct InsideNodeFeature {
    seq: String,        // k-mer sequence
    cov: usize,         // KC
    indeg: usize,
    outdeg: usize,
    is_branch: bool,
    is_articulation: bool, // undirected articulation in induced subgraph
    dist_start: usize,     // BFS hops from start
    dist_end: usize,       // BFS hops to end
    centrality_hint: f64,  // 1/(1+|ds-de|)
    cov_z: f64,            // z-score in induced
}

#[inline]
fn parse_u128(s: &str) -> u128 {
    s.parse::<u128>()
        .unwrap_or_else(|_| panic!("Expected numeric node id, got: {s}"))
}

#[inline]
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

#[inline]
fn canonical(seq: &str) -> String {
    let rc = revcomp(seq);
    if seq <= &rc { seq.to_string() } else { rc }
}

fn parse_gfa(
    path: &str,
) -> (HashMap<u128, Segment>, Vec<Link>, HashMap<(u128, u128), usize>) {
    let f = File::open(path).expect("Cannot open GFA");
    let rdr = BufReader::new(f);

    let mut segs: HashMap<u128, Segment> = HashMap::new();
    let mut links: Vec<Link> = Vec::new();
    let mut edge_cov: HashMap<(u128, u128), usize> = HashMap::new();

    for line in rdr.lines().flatten() {
        if line.is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        match cols.get(0).copied() {
            Some("S") => {
                if cols.len() < 3 { continue; }
                let id: u128 = cols[1].parse().unwrap_or_else(|_| panic!("Non-numeric segment id: {}", cols[1]));
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
                if cols.len() < 6 { continue; }
                let from: u128 = cols[1].parse().unwrap_or_else(|_| panic!("Non-numeric link from: {}", cols[1]));
                let to: u128   = cols[3].parse().unwrap_or_else(|_| panic!("Non-numeric link to: {}", cols[3]));
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
        if let Some(v) = succ.get_mut(&l.from) { v.push(l.to); }
        if let Some(v) = pred.get_mut(&l.to)   { v.push(l.from); }
    }
    (succ, pred)
}

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

    let link_cov = |u: u128, v: u128| -> usize {
        let ec = edge_cov.get(&(u, v)).copied().unwrap_or(0);
        if ec > 0 { ec } else {
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
                if !nodes.contains(&v0) { continue; }
                let mut chain = vec![u, v0];
                let mut cur = v0;

                let mut cov_sum: usize = link_cov(u, v0);
                let mut cov_min: usize = link_cov(u, v0);
                let mut edges_cnt = 1usize;

                while nodes.contains(&cur) && indeg[&cur] == 1 && outdeg[&cur] == 1 {
                    let next = succ.get(&cur).and_then(|v| v.first()).copied();
                    if let Some(nxt) = next {
                        if !nodes.contains(&nxt) { break; }
                        let ec = link_cov(cur, nxt);
                        cov_min = cmp::min(cov_min, ec);
                        cov_sum += ec;
                        edges_cnt += 1;
                        chain.push(nxt);
                        cur = nxt;
                        if indeg[&cur] != 1 || outdeg[&cur] != 1 { break; }
                    } else { break; }
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

#[inline]
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
    if chain.is_empty() { return None; }
    let mut seq = segments[&chain[0]].seq.clone();
    let mut ok = true;
    for node in chain.iter().skip(1) {
        let s = &segments[node].seq;
        if !append_next_kmer(&mut seq, s, k) { ok = false; break; }
    }
    if ok { return Some(seq); }
    let mut seq2 = revcomp(&segments[&chain[0]].seq);
    for node in chain.iter().skip(1) {
        let s = &segments[node].seq;
        if !append_next_kmer(&mut seq2, s, k) { return None; }
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
        if !append_next_kmer(assembled, s, k) { return false; }
    }
    true
}

/// Enumerate paths with a **deadline**. If deadline is hit, `timed_out=true`.
fn enumerate_paths_bidir(
    start: u128,
    end: u128,
    edges: &[SuperEdge],
    segments: &HashMap<u128, Segment>,
    k: usize,
    max_paths: usize,
    deadline: Instant,
) -> (Vec<Vec<usize>>, Vec<String>, bool, bool) {
    let from_map = from_index(edges);

    fn enumerate_from(
        s: u128,
        t: u128,
        edges: &[SuperEdge],
        from_map: &HashMap<u128, Vec<usize>>,
        segments: &HashMap<u128, Segment>,
        k: usize,
        max_paths: usize,
        deadline: Instant,
    ) -> (Vec<Vec<usize>>, Vec<String>, bool) {
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
            if Instant::now() >= deadline { return (paths, seqs, true); }
            if paths.len() >= max_paths { break; }
            if at == t {
                paths.push(idxs.clone());
                seqs.push(seq.clone());
                continue;
            }
            if idxs.len() > guard_max_steps { continue; }
            if let Some(eidxs) = from_map.get(&at) {
                for &ei in eidxs {
                    if Instant::now() >= deadline { return (paths, seqs, true); }
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
        (paths, seqs, false)
    }

    let (p_fwd, s_fwd, to1) = enumerate_from(start, end, edges, &from_map, segments, k, max_paths, deadline);
    if !p_fwd.is_empty() || to1 { return (p_fwd, s_fwd, false, to1); }
    let (p_rev, s_rev, to2) = enumerate_from(end, start, edges, &from_map, segments, k, max_paths, deadline);
    (p_rev, s_rev, true, to2)
}

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
    forward: bool,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    indeg: &HashMap<u128, usize>,
    outdeg: &HashMap<u128, usize>,
    nodes_cap: usize,
) -> Vec<u128> {
    let mut chain: Vec<u128> = vec![seed];
    let mut cur = seed;
    loop {
        if chain.len() >= nodes_cap { break; }
        let (nexts, cond) = if forward {
            (succ.get(&cur).cloned().unwrap_or_default(),
             outdeg.get(&cur).copied().unwrap_or(0) == 1)
        } else {
            (pred.get(&cur).cloned().unwrap_or_default(),
             indeg.get(&cur).copied().unwrap_or(0) == 1)
        };
        if !cond || nexts.len() != 1 { break; }
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
    if nodes.is_empty() { return Some(("".to_string(), 0)); }
    let mut seq = segments[&nodes[0]].seq.clone();
    for node in nodes.iter().skip(1) {
        let s = &segments[node].seq;
        if !append_next_kmer(&mut seq, s, k) {
            let mut seq2 = revcomp(&segments[&nodes[0]].seq);
            for node2 in nodes.iter().skip(1) {
                let s2 = &segments[node2].seq;
                if !append_next_kmer(&mut seq2, s2, k) { return None; }
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
        start, false, succ, pred, indeg, outdeg,
        max_bp / cmp::max(k, 1) + 2,
    );
    let left_nodes: Vec<u128> = left_nodes_rev.into_iter().rev().collect();
    let (left_seq_full, _) = sequence_for_nodes(&left_nodes, segments, k).unwrap_or(("".to_string(), 0));
    let right_nodes = gather_unique_chain(
        end, true, succ, pred, indeg, outdeg,
        max_bp / cmp::max(k, 1) + 2,
    );
    let (right_seq_full, _) = sequence_for_nodes(&right_nodes, segments, k).unwrap_or(("".to_string(), 0));

    let left_bp = if left_seq_full.len() > k { &left_seq_full[..left_seq_full.len() - k] } else { "" };
    let right_bp = if right_seq_full.len() > k { &right_seq_full[k..] } else { "" };
    (left_bp.to_string(), right_bp.to_string())
}

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
            .map(|(n, s)| Chrom { name: n, seq: s.clone(), rc_seq: revcomp(&s) })
            .collect();
        Self { chroms }
    }
}

#[derive(Clone, Copy, Debug)]
struct Pos { chrom: usize, pos: usize, rev: bool }

fn all_kmer_pos(reftext: &RefText, kmer: &str) -> Vec<Pos> {
    fn scan(hay: &str, needle: &str, chrom: usize, rev: bool, out: &mut Vec<Pos>) {
        if needle.is_empty() || hay.len() < needle.len() { return; }
        let nb = needle.as_bytes();
        let hb = hay.as_bytes();
        let m = nb.len();
        let mut i = 0usize;
        while i + m <= hb.len() {
            if &hb[i..i+m] == nb { out.push(Pos { chrom, pos: i, rev }); }
            i += 1;
        }
    }
    let mut out = Vec::new();
    let rc = revcomp(kmer);
    for (i, c) in reftext.chroms.iter().enumerate() {
        scan(&c.seq, kmer, i, false, &mut out);
        scan(&c.rc_seq, &rc, i, true, &mut out);
    }
    out
}

fn ref_substr_between(reftext: &RefText, a: Pos, b: Pos, k: usize) -> Option<String> {
    if a.chrom != b.chrom || a.rev != b.rev { return None; }
    let c = &reftext.chroms[a.chrom];
    if !a.rev {
        if a.pos > b.pos { return None; }
        let start = a.pos;
        let end_excl = b.pos + k;
        if end_excl <= c.seq.len() { return Some(c.seq[start..end_excl].to_string()); }
    } else {
        if a.pos > b.pos { return None; }
        let start = a.pos;
        let end_excl = b.pos + k;
        if end_excl <= c.rc_seq.len() { return Some(c.rc_seq[start..end_excl].to_string()); }
    }
    None
}

fn minimal_unique_extension(
    reftext: &RefText,
    base_kmer: &str,
    extend_right: bool,
    max_ext: usize,
) -> Option<(usize, Pos)> {
    let mut occ = all_kmer_pos(reftext, base_kmer);
    if occ.is_empty() { return None; }
    if occ.len() == 1 { return Some((0, occ[0])); }

    for ext in 1..=max_ext {
        use std::collections::HashMap;
        let mut groups: HashMap<(usize,u8,bool), Vec<Pos>> = HashMap::new();
        for p in occ.iter().copied() {
            let c = &reftext.chroms[p.chrom];
            let next_base = if !p.rev {
                if extend_right {
                    c.seq.as_bytes().get(p.pos + base_kmer.len() + (ext-1)).copied()
                } else {
                    if p.pos >= ext { c.seq.as_bytes().get(p.pos - ext).copied() } else { None }
                }
            } else {
                if extend_right {
                    c.rc_seq.as_bytes().get(p.pos + base_kmer.len() + (ext-1)).copied()
                } else {
                    if p.pos >= ext { c.rc_seq.as_bytes().get(p.pos - ext).copied() } else { None }
                }
            };
            if let Some(b) = next_base {
                groups.entry((p.chrom,b,p.rev)).or_default().push(p);
            }
        }
        let mut unique_groups: Vec<Pos> = Vec::new();
        for v in groups.values() {
            if v.len() == 1 { unique_groups.push(v[0]); }
        }
        if unique_groups.len() == 1 {
            return Some((ext, unique_groups[0]));
        }
        occ = groups.values().flat_map(|v| v.clone()).collect();
        if occ.len() == 1 { return Some((ext, occ[0])); }
    }
    None
}

fn label_with_ref_walk(
    candidate_seqs: &[String],
    start_seq: &str,
    end_seq: &str,
    reftext: &RefText,
    k: usize,
    try_mue: bool,
    deadline: Instant,
) -> (Option<usize>, String) {
    let mut starts = all_kmer_pos(reftext, start_seq);
    let mut ends   = all_kmer_pos(reftext, end_seq);

    if Instant::now() >= deadline { return (None, "timeout".into()); }

    if try_mue && (starts.len() > 2 || ends.len() > 2) {
        if let Some((_e, p)) = minimal_unique_extension(reftext, start_seq, true, 4096) {
            starts = vec![p];
        }
        if Instant::now() >= deadline { return (None, "timeout".into()); }
        if let Some((_e, p)) = minimal_unique_extension(reftext, end_seq, false, 4096) {
            ends = vec![p];
        }
    }

    let mut matches: Vec<usize> = Vec::new();
    for s in &starts {
        if Instant::now() >= deadline { return (None, "timeout".into()); }
        for e in &ends {
            if Instant::now() >= deadline { return (None, "timeout".into()); }
            if s.chrom != e.chrom || s.rev != e.rev { continue; }
            if let Some(ref_sub) = ref_substr_between(reftext, *s, *e, k) {
                let refb = ref_sub.as_bytes();
                for (i, cand) in candidate_seqs.iter().enumerate() {
                    if Instant::now() >= deadline { return (None, "timeout".into()); }
                    let cb = cand.as_bytes();
                    if refb.len() == cb.len() && refb == cb {
                        matches.push(i);
                    } else {
                        let rc = revcomp(cand);
                        if refb.len() == rc.len() && refb == rc.as_bytes() {
                            matches.push(i);
                        }
                    }
                }
            }
        }
    }
    matches.sort_unstable();
    matches.dedup();

    if matches.len() == 1 {
        (Some(matches[0]), "unique_ref_walk".into())
    } else if matches.is_empty() {
        (None, "no_ref_match".into())
    } else {
        (None, "ambiguous_ref_multi_hit".into())
    }
}

fn label_with_adaptive_flanks_lenient_ac(
    candidate_seqs: &[String],
    start: u128,
    end: u128,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    segments: &HashMap<u128, Segment>,
    k: usize,
    reftext: &RefText,
    deadline: Instant,
) -> (Option<usize>, usize, bool, bool) {
    use std::collections::{HashMap, HashSet};

    // RC-class collapse
    let mut class_map: HashMap<String, (usize, Vec<usize>)> = HashMap::new();
    for (i, s) in candidate_seqs.iter().enumerate() {
        let key = canonical(s);
        class_map.entry(key).and_modify(|(_, v)| v.push(i)).or_insert((i, vec![i]));
    }
    let classes: Vec<(String, usize, Vec<usize>)> =
        class_map.into_iter().map(|(k,(rep,mem))| (k,rep,mem)).collect();

    if Instant::now() >= deadline { return (None, 0, false, true); }

    let max_flank_bp = 4_096usize;
    let (indeg, outdeg) = compute_global_deg(succ, pred);
    let (left_full, right_full) = compute_flanks(start, end, succ, pred, &indeg, &outdeg, segments, k, max_flank_bp);
    let flanks: &[usize] = &[0, 32, 64, 128, 256, 512, 1024, 2048, 4096];

    for d in flanks {
        if Instant::now() >= deadline { return (None, d * 2, false, true); }
        let left_tail  = if left_full.len()  > *d { &left_full[left_full.len() - *d..] } else { &left_full[..] };
        let right_head = if right_full.len() > *d { &right_full[..*d] } else { &right_full[..] };

        let mut patterns: Vec<String> = Vec::with_capacity(classes.len() * 2);
        for (class_seq, _, _) in classes.iter() {
            patterns.push(format!("{left}{mid}{right}", left = left_tail, mid = class_seq, right = right_head));
            patterns.push(format!("{left}{mid}{right}", left = left_tail, mid = revcomp(class_seq), right = right_head));
        }

        if Instant::now() >= deadline { return (None, d * 2, false, true); }

        let ac = AhoCorasickBuilder::new()
            .match_kind(MatchKind::LeftmostFirst)
            .build(patterns.iter().map(|s| s.as_str()))
            .expect("Failed to build Aho–Corasick automaton");

        let mut class_hits: Vec<usize> = vec![0; classes.len()];
        let mut seen: HashSet<(usize, usize, usize)> = HashSet::new();

        for (cx, chrom) in reftext.chroms.iter().enumerate() {
            if Instant::now() >= deadline { return (None, d * 2, false, true); }
            for m in ac.find_iter(&chrom.seq) {
                let pid = m.pattern().as_usize();
                let class_ix = pid / 2;
                let dir_mod  = pid % 2;
                let key = (cx, m.start(), dir_mod);
                if seen.insert(key) {
                    class_hits[class_ix] += 1;
                }
            }
        }
        let winners: Vec<usize> = class_hits.iter().enumerate().filter_map(|(i,&h)| (h>0).then_some(i)).collect();
        if winners.len() == 1 {
            let ci = winners[0];
            let rep_idx = classes[ci].1;
            let topology_tie = classes[ci].2.len() > 1;
            return (Some(rep_idx), d * 2, topology_tie, false);
        }
    }
    (None, max_flank_bp * 2, false, false)
}

fn count_bubbles_file(path: &str) -> usize {
    let in_file = File::open(path).expect("open bubbles json/jsonl");
    let reader = BufReader::new(in_file);
    let stream = Deserializer::from_reader(reader).into_iter::<Value>();
    let mut total = 0usize;
    for val in stream {
        let Ok(v) = val else { continue };
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

// ---------- Strategy A helpers (unsupervised) ----------

fn induced_degrees(
    nodes: &HashSet<u128>,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
) -> (HashMap<u128, usize>, HashMap<u128, usize>) {
    let mut indeg = HashMap::new();
    let mut outdeg = HashMap::new();
    for &u in nodes {
        let o = succ.get(&u).map(|v| v.iter().filter(|x| nodes.contains(x)).count()).unwrap_or(0);
        let i = pred.get(&u).map(|v| v.iter().filter(|x| nodes.contains(x)).count()).unwrap_or(0);
        indeg.insert(u, i);
        outdeg.insert(u, o);
    }
    (indeg, outdeg)
}

fn bfs_distances_from(
    seeds: &[u128],
    nodes: &HashSet<u128>,
    adj: &HashMap<u128, Vec<u128>>,
) -> HashMap<u128, usize> {
    let mut dist: HashMap<u128, usize> = HashMap::new();
    let mut dq = VecDeque::new();
    for &s in seeds {
        if nodes.contains(&s) {
            dist.insert(s, 0);
            dq.push_back(s);
        }
    }
    while let Some(u) = dq.pop_front() {
        if let Some(vs) = adj.get(&u) {
            for &v in vs {
                if !nodes.contains(&v) { continue; }
                if dist.contains_key(&v) { continue; }
                dist.insert(v, dist[&u] + 1);
                dq.push_back(v);
            }
        }
    }
    dist
}

fn articulation_points_undirected(
    nodes: &HashSet<u128>,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
) -> HashSet<u128> {
    // build undirected adjacency
    let mut und: HashMap<u128, Vec<u128>> = HashMap::new();
    for &u in nodes {
        und.entry(u).or_default();
        for &v in succ.get(&u).unwrap_or(&Vec::new()) {
            if nodes.contains(&v) {
                und.entry(u).or_default().push(v);
                und.entry(v).or_default().push(u);
            }
        }
        for &p in pred.get(&u).unwrap_or(&Vec::new()) {
            if nodes.contains(&p) {
                und.entry(u).or_default().push(p);
                und.entry(p).or_default().push(u);
            }
        }
    }
    let mut time = 0usize;
    let mut disc: HashMap<u128, usize> = HashMap::new();
    let mut low: HashMap<u128, usize> = HashMap::new();
    let mut parent: HashMap<u128, Option<u128>> = HashMap::new();
    let mut ap: HashSet<u128> = HashSet::new();

    fn dfs(
        u: u128,
        und: &HashMap<u128, Vec<u128>>,
        time: &mut usize,
        disc: &mut HashMap<u128, usize>,
        low: &mut HashMap<u128, usize>,
        parent: &mut HashMap<u128, Option<u128>>,
        ap: &mut HashSet<u128>,
    ) {
        *time += 1;
        let t = *time;
        disc.insert(u, t);
        low.insert(u, t);
        let mut children = 0usize;

        if let Some(nei) = und.get(&u) {
            for &v in nei {
                if !disc.contains_key(&v) {
                    children += 1;
                    parent.insert(v, Some(u));
                    dfs(v, und, time, disc, low, parent, ap);
                    let low_u = *low.get(&u).unwrap();
                    let low_v = *low.get(&v).unwrap();
                    low.insert(u, low_u.min(low_v));
                    let p = parent.get(&u).copied().flatten();
                    if p.is_none() && children > 1 { ap.insert(u); }
                    if p.is_some() && low_v >= *disc.get(&u).unwrap() { ap.insert(u); }
                } else if Some(v) != parent.get(&u).copied().flatten() {
                    let low_u = *low.get(&u).unwrap();
                    let disc_v = *disc.get(&v).unwrap();
                    low.insert(u, low_u.min(disc_v));
                }
            }
        }
    }

    for &u in und.keys() {
        if !disc.contains_key(&u) {
            parent.insert(u, None);
            dfs(u, &und, &mut time, &mut disc, &mut low, &mut parent, &mut ap);
        }
    }
    ap
}

fn mean_std(vals: &[f64]) -> (f64, f64) {
    if vals.is_empty() { return (0.0, 1.0); }
    let m = vals.iter().sum::<f64>() / vals.len() as f64;
    let v = vals.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / vals.len() as f64;
    (m, v.sqrt().max(1e-9))
}

fn select_inside_nodes_meta(
    start: u128,
    end: u128,
    original_inside: &HashSet<u128>,
    induced_nodes: &HashSet<u128>,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    segments: &HashMap<u128, Segment>,
    top_k: usize,
) -> Vec<InsideNodeFeature> {
    let mut cand: Vec<u128> = original_inside
        .iter().copied()
        .filter(|&u| induced_nodes.contains(&u) && u != start && u != end)
        .collect();
    if cand.is_empty() { return Vec::new(); }

    let (indeg, outdeg) = induced_degrees(induced_nodes, succ, pred);
    let ds = bfs_distances_from(&[start], induced_nodes, succ);
    let de = bfs_distances_from(&[end], induced_nodes, pred);
    let aps = articulation_points_undirected(induced_nodes, succ, pred);

    let covs: Vec<f64> = induced_nodes
        .iter()
        .map(|u| segments.get(u).map(|s| s.cov as f64).unwrap_or(0.0))
        .collect();
    let (cov_mu, cov_sigma) = mean_std(&covs);

    #[derive(Clone)]
    struct Row { u: u128, score: f64, feat: InsideNodeFeature }

    let mut rows: Vec<Row> = Vec::with_capacity(cand.len());
    for u in cand.drain(..) {
        let cov = segments.get(&u).map(|s| s.cov).unwrap_or(0);
        let indeg_u = *indeg.get(&u).unwrap_or(&0);
        let outdeg_u = *outdeg.get(&u).unwrap_or(&0);
        let is_branch = indeg_u > 1 || outdeg_u > 1;
        let is_ap = aps.contains(&u);

        let d_s = *ds.get(&u).unwrap_or(&usize::MAX);
        let d_e = *de.get(&u).unwrap_or(&usize::MAX);
        let centrality_hint = if d_s == usize::MAX || d_e == usize::MAX {
            0.0
        } else {
            1.0 / (1.0 + ((d_s as isize - d_e as isize).abs() as f64))
        };

        let cov_z = if cov_sigma > 0.0 { (cov as f64 - cov_mu) / cov_sigma } else { 0.0 };
        let cov_contrast = cov_z.abs();

        let score =
            2.5 * (is_branch as i32 as f64) +
            2.0 * (is_ap as i32 as f64) +
            1.2 * centrality_hint +
            0.8 * cov_contrast;

        let seq = segments.get(&u).map(|s| s.seq.clone()).unwrap_or_default();

        rows.push(Row {
            u,
            score,
            feat: InsideNodeFeature {
                seq,
                cov,
                indeg: indeg_u,
                outdeg: outdeg_u,
                is_branch,
                is_articulation: is_ap,
                dist_start: if d_s == usize::MAX { 0 } else { d_s },
                dist_end:   if d_e == usize::MAX { 0 } else { d_e },
                centrality_hint,
                cov_z,
            }
        });
    }

    rows.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score).unwrap()
            .then_with(|| (b.feat.indeg + b.feat.outdeg).cmp(&(a.feat.indeg + a.feat.outdeg)))
            .then_with(|| b.feat.cov.cmp(&a.feat.cov))
            .then_with(|| a.u.cmp(&b.u))
    });

    let mut out: Vec<InsideNodeFeature> = Vec::new();
    let mut per_seq: HashMap<String, usize> = HashMap::new();
    for r in rows {
        if out.len() >= top_k { break; }
        let c = per_seq.entry(r.feat.seq.clone()).or_insert(0);
        if *c >= 2 { continue; }
        *c += 1;
        out.push(r.feat);
    }
    out
}

// -------------- build_record_json with Strategy A inside nodes --------------

fn build_record_json(
    bubble: &Bubble,
    k: usize,
    succ: &HashMap<u128, Vec<u128>>,
    pred: &HashMap<u128, Vec<u128>>,
    segments: &HashMap<u128, Segment>,
    edge_cov: &HashMap<(u128, u128), usize>,
    reftext: &RefText,
    radius_context: usize,
    flank_kmers: usize,
    max_paths: usize,
    deadline: Instant,
) -> Option<String> {
    if bubble.ends.len() != 2 { return None; }
    let start: u128 = parse_u128(&bubble.ends[0]);
    let end: u128   = parse_u128(&bubble.ends[1]);

    let inside: HashSet<u128> = bubble.inside.iter().map(|s| parse_u128(s)).collect();
    let nodes: HashSet<u128> =
        collect_nodes_for_bubble_context(&inside, start, end, succ, pred, radius_context);

    if Instant::now() >= deadline {
        return Some(serde_json::to_string(&BubbleRecord {
            bubble_id: bubble.id, k,
            start_seq: segments.get(&start).map(|s| s.seq.clone()).unwrap_or_default(),
            end_seq: segments.get(&end).map(|s| s.seq.clone()).unwrap_or_default(),
            upstream: vec![], downstream: vec![],
            upstream_nodes: vec![], downstream_nodes: vec![],
            nodes: vec![], edges: vec![], inside_nodes: vec![], paths: vec![],
            label_path: None, label_reason: Some("timeout".into()),
        }).unwrap());
    }

    let (superedges0, supernodes0) = compact_with_chains(&nodes, succ, pred, edge_cov, segments);
    let mut superedges = split_edges_at(start, &superedges0);
    superedges = split_edges_at(end, &superedges);

    let (paths_vec, seqs_vec, _reversed, to_enum) =
        enumerate_paths_bidir(start, end, &superedges, segments, k, max_paths, deadline);

    let upstream_ids: Vec<u128> = collect_flanks_greedy(
        start, false, succ, pred, segments, edge_cov, flank_kmers
    );
    let downstream_ids: Vec<u128> = collect_flanks_greedy(
        end, true, succ, pred, segments, edge_cov, flank_kmers
    );

    let top_inside_k: usize = 16;
    let inside_nodes = select_inside_nodes_meta(
        start, end, &inside, &nodes, succ, pred, segments, top_inside_k
    );

    let mut rec = BubbleRecord {
        bubble_id: bubble.id,
        k,
        start_seq: segments.get(&start).map(|s| s.seq.clone()).unwrap_or_default(),
        end_seq:   segments.get(&end).map(|s| s.seq.clone()).unwrap_or_default(),
        upstream: upstream_ids
            .iter().filter_map(|nid| segments.get(nid).map(|s| s.seq.clone())).collect(),
        downstream: downstream_ids
            .iter().filter_map(|nid| segments.get(nid).map(|s| s.seq.clone())).collect(),
        upstream_nodes: upstream_ids
            .iter().filter_map(|nid| segments.get(nid))
            .map(|s| NodeFeature { seq: s.seq.clone(), cov: s.cov }).collect(),
        downstream_nodes: downstream_ids
            .iter().filter_map(|nid| segments.get(nid))
            .map(|s| NodeFeature { seq: s.seq.clone(), cov: s.cov }).collect(),
        nodes: Vec::new(),
        edges: Vec::new(),
        inside_nodes,
        paths: paths_vec.clone(),
        label_path: None,
        label_reason: None,
    };

    for &nid in &supernodes0 {
        if let Some(seg) = segments.get(&nid) {
            rec.nodes.push(NodeFeature { seq: seg.seq.clone(), cov: seg.cov });
        }
    }
    for e in superedges.iter() {
        let len_nodes = e.chain.len();
        let len_bp = k + len_nodes.saturating_sub(1);
        let source_seq = segments.get(&e.from).map(|s| s.seq.clone()).unwrap_or_default();
        let target_seq = segments.get(&e.to).map(|s| s.seq.clone()).unwrap_or_default();
        rec.edges.push(EdgeFeature {
            source_seq, target_seq, len_nodes, len_bp,
            cov_min: e.cov_min, cov_mean: e.cov_mean,
        });
    }

    if to_enum || Instant::now() >= deadline {
        rec.label_path = None;
        rec.label_reason = Some("timeout".into());
        return Some(serde_json::to_string(&rec).unwrap());
    }

    if !seqs_vec.is_empty() {
        let (lab1, why1) = label_with_ref_walk(
            &seqs_vec, &segments[&start].seq, &segments[&end].seq,
            reftext, k, true, deadline,
        );
        if why1 == "timeout" {
            rec.label_reason = Some("timeout".into());
            return Some(serde_json::to_string(&rec).unwrap());
        }
        if let Some(ix) = lab1 {
            rec.label_path = Some(ix);
            rec.label_reason = Some(why1);
        } else {
            let (lab2, _bp, seq_tie, timed_out) = label_with_adaptive_flanks_lenient_ac(
                &seqs_vec, start, end, succ, pred, segments, k, reftext, deadline
            );
            if timed_out {
                rec.label_reason = Some("timeout".into());
            } else if lab2.is_some() && !seq_tie {
                rec.label_path = lab2;
                rec.label_reason = Some("adaptive_kmp".into());
            } else if lab2.is_some() && seq_tie {
                rec.label_path = None;
                rec.label_reason = Some("adaptive_kmp_sequence_tie".into());
            } else {
                rec.label_reason = Some(why1);
            }
        }
    } else {
        rec.label_reason = Some("no_paths".into());
    }

    Some(serde_json::to_string(&rec).unwrap())
}
fn collect_flanks_greedy(
    seed: u128,
    forward: bool,
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
            if seen.contains(&nx) { break; } // avoid tiny cycles
            out.push(nx);
            seen.insert(nx);
            cur = nx;
        } else { break; }
    }
    out
}

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
    if cand.is_empty() { return None; }
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
        if d == radius { continue; }
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
        if d == radius { continue; }
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

    // Tunables (safe defaults; feel free to raise for more recall)
    let radius_context: usize = 12;   // nodes to include around bubble for compaction
    let flank_kmers: usize = 5;       // preset number of upstream/downstream kmers to export
    let max_paths: usize = 2048;      // candidate cap per bubble
    const BATCH: usize = 200;         // bubbles per parallel batch (tweak for memory/CPU)

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
        let v = match val { Ok(v) => v, Err(_) => continue };

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
                for bc in bcs { bubbles.extend(bc.bubbles); }
            } else if let Ok(bs) = serde_json::from_value::<Vec<Bubble>>(v.clone()) {
                bubbles.extend(bs);
            } else {
                panic!("Unsupported array JSON structure for bubbles");
            }
        } else if let Ok(bc) = serde_json::from_value::<BubbleChain>(v.clone()) {
            bubbles.extend(bc.bubbles);
        } else {
            continue;
        }

        for b in bubbles {
            batch.push(b);
            if batch.len() == BATCH {
                // Move out the batch to avoid borrow issues, then process in parallel.
                let curr_batch: Vec<Bubble> = std::mem::take(&mut batch);
                let pb_for_threads = pb.clone();
                let mut lines: Vec<(usize, String)> = curr_batch
                    .par_iter()
                    .enumerate()
                    .map_init(|| pb_for_threads.clone(), |pbl, (i, bub)| {
                        let deadline = Instant::now() + Duration::from_secs(PER_BUBBLE_TIMEOUT_SECS);
                        let out = build_record_json(
                            bub, k, &succ, &pred, &segments, &edge_cov, &reftext,
                            radius_context, flank_kmers, max_paths, deadline,
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
                total_written += curr_batch.len();
            }
        }
    }

    // flush tail
    if !batch.is_empty() {
        let curr_batch: Vec<Bubble> = std::mem::take(&mut batch);
        let pb_for_threads = pb.clone();
        let mut lines: Vec<(usize, String)> = curr_batch
            .par_iter()
            .enumerate()
            .map_init(|| pb_for_threads.clone(), |pbl, (i, bub)| {
                let deadline = Instant::now() + Duration::from_secs(PER_BUBBLE_TIMEOUT_SECS);
                let out = build_record_json(
                    bub, k, &succ, &pred, &segments, &edge_cov, &reftext,
                    radius_context, flank_kmers, max_paths, deadline,
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
        total_written += curr_batch.len();
    }

    pb.finish_with_message(format!("Processed {} bubbles", total_bubbles));
    eprintln!("Done: wrote {} hyperbubble records", total_written);
}
