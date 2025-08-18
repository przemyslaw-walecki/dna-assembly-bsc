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
//! - **Fallback labeler:** adaptive flank extension + exact **KMP** search on both strands.
//! - **ID-free output:** upstream/downstream and endpoints are sequences.
//! - Coverage: `cov_min` / `cov_mean` prefer **EC:i** per link; fallback to `min(KC:u, KC:v)`.
//! - **Parallelized** with Rayon (batched to bound memory); global progress via `indicatif`.
//!
//! ## Usage
//! ```text
//! cargo run --release -- <graph.gfa> <bubbles.json|jsonl> <reference.fasta> <out.jsonl>
//! ```
//!
//! ## Output (JSONL per bubble)
//! Each line is a `BubbleRecord`:
//! ```json
//! {
//!   "bubble_id": 123,
//!   "k": 31,
//!   "start_seq": "ATGC...",
//!   "end_seq":   "TGCA...",
//!   "upstream": ["..."], "downstream": ["..."],
//!   "upstream_nodes": [{"seq":"...","cov":42}], "downstream_nodes": [{"seq":"...","cov":17}],
//!   "nodes": [{"seq":"...","cov":...}],
//!   "edges": [
//!     {"source_seq":"...","target_seq":"...","len_nodes":7,"len_bp":37,"cov_min":10,"cov_mean":12.5}
//!   ],
//!   "paths": [[0,2,5], [1,4]],
//!   "label_path": 0,
//!   "label_reason": "unique_ref_walk" // or "adaptive_kmp", "ambiguous_ref_multi_hit", "no_ref_match"
//! }
//! ```
//!
//! ## Notes
//! - `k` is inferred from the first segment’s sequence length in the GFA.
//! - JSON input supports `BubbleChain` collections or raw `Bubble` arrays (JSON or JSONL).
//! - Batching keeps memory usage stable during parallel processing.

use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::{Deserializer, Value};
use std::cmp;
use std::collections::{HashMap, HashSet, VecDeque};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/// Segment (GFA `S` record).
///
/// - `id`: numeric segment identifier (required to be numeric in this tool).
/// - `seq`: segment sequence (uppercase A/C/G/T/N).
/// - `cov`: node coverage (from `KC:i` if present; otherwise 0).
#[derive(Debug)]
struct Segment {
    id: u128,
    seq: String,
    cov: usize, // KC
}

/// Link (GFA `L` record).
///
/// - `from` → `to`: numeric segment ids.
/// - `cov`: edge/link coverage (from `EC:i` if present; otherwise 0).
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
///
/// - `id`: stable bubble identifier from upstream discovery.
/// - `ends`: two endpoint node identifiers as strings (numeric content expected).
/// - `inside`: node identifiers (strings) forming the internal bubble subgraph.
#[derive(Deserialize)]
struct Bubble {
    id: usize,
    ends: Vec<String>,   // [start,end]
    inside: Vec<String>, // node ids
}

/// A compacted path (“superedge”/unitig) across nodes with coverage aggregates.
///
/// `chain` includes endpoints `from` and `to`. The `chain` is **not serialized**;
/// downstream JSON only exposes endpoint sequences.
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
    /// Bubble identifier (from input).
    bubble_id: usize,
    /// *k*-mer size inferred from GFA segment sequence length.
    k: usize,

    /// Endpoint sequences.
    start_seq: String,
    end_seq: String,

    /// Greedily selected upstream/downstream neighbor sequences (preset count).
    upstream: Vec<String>,
    downstream: Vec<String>,

    /// Upstream/downstream neighbor nodes with coverage.
    upstream_nodes: Vec<NodeFeature>,
    downstream_nodes: Vec<NodeFeature>,

    /// All supernodes’ sequences (context) with coverage.
    nodes: Vec<NodeFeature>,

    /// Superedges serialized by endpoint sequences and coverage aggregates.
    edges: Vec<EdgeFeature>,

    /// Candidate paths as lists of edge indices into `edges`.
    paths: Vec<Vec<usize>>,
    /// Chosen label path index (if any).
    label_path: Option<usize>,
    /// Labeling rationale.
    /// One of: `"unique_ref_walk"`, `"adaptive_kmp"`, `"ambiguous_ref_multi_hit"`, `"no_ref_match"`.
    label_reason: Option<String>,
}

/// Node-level feature (sequence + coverage).
#[derive(Serialize)]
struct NodeFeature {
    seq: String,
    cov: usize,
}

/// Edge-level feature between two endpoint sequences with coverage aggregates.
#[derive(Serialize)]
struct EdgeFeature {
    source_seq: String,
    target_seq: String,
    /// Number of nodes in this superedge path (`chain.len()`).
    len_nodes: usize,
    /// Total base pairs spanned assuming overlapping *k* (`k + len_nodes - 1`, saturating at `k`).
    len_bp: usize,
    /// Minimum edge coverage across the chain (prefers `EC:i`, else `min(KC)`).
    cov_min: usize,
    /// Mean edge coverage across the chain.
    cov_mean: f64,
}

/// Parse a string into `u128` with a clear panic on failure.
///
/// Used to normalize string node ids from JSON into numeric form.
fn parse_u128(s: &str) -> u128 {
    s.parse::<u128>()
        .unwrap_or_else(|_| panic!("Expected numeric node id, got: {s}"))
}

/// Reverse-complement of a DNA string (A↔T, C↔G; others→N).
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

/// Parse a GFA file to collect segments, links, and a map of edge coverages.
///
/// Coverage preference:
/// - Edge coverage `EC:i` (per `L` record) if present.
/// - Otherwise fall back to `min(KC:u, KC:v)` when aggregating.
///
/// Returns:
/// 1) map of `id → Segment`
/// 2) list of `Link`
/// 3) map of `(from, to) → EC`
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
                if cols.len() < 3 { continue; }
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
                if cols.len() < 6 { continue; }
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

/// Build adjacency lists for successors and predecessors from segments and links.
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

/// Compact simple chains into `SuperEdge`s and compute coverage aggregates.
///
/// Compaction respects an induced node set `nodes`. Coverage preference is
/// `EC:i` when available; otherwise `min(KC:u, KC:v)` is used.
///
/// Returns:
/// - vector of `SuperEdge`s,
/// - list of “supernodes” (nodes with `in≠1` or `out≠1`) forming endpoints.
fn compact_with_chains(
    nodes: &HashSet<u128>,
    //
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

/// Return the position of a node in a chain (if present).
fn index_in_chain(chain: &[u128], node: u128) -> Option<usize> {
    chain.iter().position(|&x| x == node)
}

/// Split superedges at a given node, preserving coverage aggregates.
///
/// If the split node is interior, two new edges are emitted. Edges that only
/// touch at endpoints are passed through unchanged.
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

/// Build an index from `from` node → list of edge indices.
fn from_index(edges: &[SuperEdge]) -> HashMap<u128, Vec<usize>> {
    let mut m: HashMap<u128, Vec<usize>> = HashMap::new();
    for (i, e) in edges.iter().enumerate() {
        m.entry(e.from).or_default().push(i);
    }
    m
}

/// Append the *k*-mer from `next_seq` to `assembled` if it overlaps by `k-1`.
///
/// Attempts forward first; if mismatch, tries reverse-complement of `next_seq`.
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

/// Assemble a sequence for a chain of nodes (orientation-safe).
///
/// Tries starting from the first node as-is; if inconsistent, retries with its RC.
fn assemble_chain_sequence(
    chain: &[u128],
    segments: &HashMap<u128, Segment>,
    k: usize,
) -> Option<String> {
    if chain.is_empty() { return None; }
    // try direct
    let mut seq = segments[&chain[0]].seq.clone();
    let mut ok = true;
    for node in chain.iter().skip(1) {
        let s = &segments[node].seq;
        if !append_next_kmer(&mut seq, s, k) { ok = false; break; }
    }
    if ok { return Some(seq); }
    // retry with RC of first
    let mut seq2 = revcomp(&segments[&chain[0]].seq);
    for node in chain.iter().skip(1) {
        let s = &segments[node].seq;
        if !append_next_kmer(&mut seq2, s, k) { return None; }
    }
    Some(seq2)
}

/// Append a chain’s sequence to an already assembled string (orientation-safe).
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

/// Enumerate paths from `start` to `end` (and reverse if forward is empty),
/// assembling orientation-safe path sequences.
///
/// Returns `(paths, sequences, reversed)` where:
/// - `paths` are edge-index lists,
/// - `sequences` are assembled strings per path,
/// - `reversed` indicates whether reverse enumeration was used.
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
            if paths.len() >= max_paths { break; }
            if at == t {
                paths.push(idxs.clone());
                seqs.push(seq.clone());
                continue;
            }
            if idxs.len() > guard_max_steps { continue; }
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

/// Compute global in/out-degrees for the whole graph (used in labeling).
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

/// Gather a unique-degree chain from a seed by repeatedly following the sole
/// successor (or predecessor), capped by `nodes_cap`.
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

/// Assemble sequence for an ordered list of node ids. Returns assembled string
/// and extra bp beyond the first *k* (or zero on failure).
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
            // retry with RC of first
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

/// Compute left/right flanking sequences by walking unique-degree chains from
/// the endpoints up to `max_bp` (converted to node count by `k`).
///
/// Returns `(left_bp, right_bp)` excluding the overlapping *k*.
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

/// Build KMP prefix function (LPS array) for a pattern.
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

/// Find all occurrences of `pat` in `text` using KMP. Returns start indices.
fn kmp_search_all(pat: &str, text: &str) -> Vec<usize> {
    if pat.is_empty() { return vec![]; }
    let lps = kmp_build(pat);
    let pb = pat.as_bytes();
    let tb = text.as_bytes();
    let mut res = Vec::new();
    let mut i = 0usize;
    let mut j = 0usize;
    while i < tb.len() {
        if tb[i] == pb[j] {
            i += 1; j += 1;
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

/// A reference chromosome with both forward and reverse-complement sequences.
#[derive(Clone)]
struct Chrom {
    name: String,
    seq: String,
    rc_seq: String,
}

/// Reference text container for all chromosomes (both strands).
struct RefText {
    chroms: Vec<Chrom>,
}

impl RefText {
    /// Load a FASTA file into `RefText`, computing reverse-complements.
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

/// Reference position of a *k*-mer.
#[derive(Clone, Copy, Debug)]
struct Pos { chrom: usize, pos: usize, rev: bool } // rev=false => forward strand

/// Find all occurrences of a *k*-mer on both strands of the reference.
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

/// Extract the reference substring between two *k*-mer positions on the same
/// chromosome and strand, inclusive of both end *k*-mers.
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

/// Minimal Unique Extension (MUE) of a base *k*-mer on the reference.
///
/// Iteratively extends to the right (`extend_right=true`) or left (`false`)
/// until a unique occurrence is found or `max_ext` is reached. Returns the
/// extension length and the unique position if successful.
fn minimal_unique_extension(
    reftext: &RefText,
    base_kmer: &str,
    extend_right: bool, // true for start(+), false for end(+)
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

/// Primary labeler: compare candidate path sequences against the exact
/// reference substring between endpoint *k*-mers (both strands).
///
/// If endpoints are repetitive, optionally applies MUE to disambiguate.
///
/// Returns `(Some(index), "unique_ref_walk")` on a unique match; otherwise
/// `(None, reason)` where `reason` is `"no_ref_match"` or `"ambiguous_ref_multi_hit"`.
fn label_with_ref_walk(
    candidate_seqs: &[String],
    start_seq: &str,
    end_seq: &str,
    reftext: &RefText,
    k: usize,
    try_mue: bool,
) -> (Option<usize>, String) {
    let mut starts = all_kmer_pos(reftext, start_seq);
    let mut ends   = all_kmer_pos(reftext, end_seq);

    if try_mue && (starts.len() > 8 || ends.len() > 8) {
        if let Some((_e, p)) = minimal_unique_extension(reftext, start_seq, true, 4096) {
            starts = vec![p];
        }
        if let Some((_e, p)) = minimal_unique_extension(reftext, end_seq, false, 4096) {
            ends = vec![p];
        }
    }

    let mut matches: Vec<usize> = Vec::new();
    for s in &starts {
        for e in &ends {
            if s.chrom != e.chrom || s.rev != e.rev { continue; }
            if let Some(ref_sub) = ref_substr_between(reftext, *s, *e, k) {
                let refb = ref_sub.as_bytes();
                for (i, cand) in candidate_seqs.iter().enumerate() {
                    let cb = cand.as_bytes();
                    if refb.len() == cb.len() && refb == cb {
                        matches.push(i);
                    } else {
                        // also try reverse-complement of candidate
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

/// Fallback labeler: adaptively extend flanks and search with exact KMP on
/// both strands. Chooses the unique path that yields exactly one global hit.
///
/// Returns `(Some(index), used_flank_bp)` on success; otherwise `(None, max_flank_bp*2)`.
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
        } else { &left_full[..] };
        let right_head = if right_full.len() > delta {
            &right_full[..delta]
        } else { &right_full[..] };

        let mut hits_per_path: Vec<usize> = Vec::with_capacity(candidate_seqs.len());
        for seq in candidate_seqs {
            // try both strands by checking seq and revcomp(seq) as separate queries
            let extended_f = format!("{left}{mid}{right}", left = left_tail, mid = seq, right = right_head);
            let extended_r = format!("{left}{mid}{right}", left = left_tail, mid = revcomp(seq), right = right_head);
            let mut total = 0usize;
            for chrom in reftext.chroms.iter() {
                total += kmp_search_all(&extended_f, &chrom.seq).len();
                total += kmp_search_all(&extended_f, &chrom.rc_seq).len();
                total += kmp_search_all(&extended_r, &chrom.seq).len();
                total += kmp_search_all(&extended_r, &chrom.rc_seq).len();
            }
            hits_per_path.push(total);
        }

        let mut unique_idx: Option<usize> = None;
        let mut bad = false;
        for (i, &h) in hits_per_path.iter().enumerate() {
            match h {
                0 => {}
                1 => { if unique_idx.is_none() { unique_idx = Some(i); } else { bad = true; break; } }
                _ => { bad = true; break; }
            }
        }
        if !bad {
            if let Some(i) = unique_idx {
                return (Some(i), delta * 2);
            }
        }
    }
    (None, max_flank_bp * 2)
}

/// Choose the “best” neighbor by (EC desc, KC desc, id asc) for determinism.
///
/// Direction controlled by `forward`: `true` uses successors, `false` predecessors.
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

/// Collect exactly `n` neighbors greedily from `seed` in the given direction,
/// avoiding tiny cycles. Used for exporting preset flanks.
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

/// Collect an induced node set for a bubble context by expanding:
/// - upstream from `start` by `radius`,
/// - downstream from `end` by `radius`,
/// and uniting with the bubble’s internal nodes.
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

/// Count total bubbles in a JSON/JSONL file without loading everything into memory.
/// Supports `BubbleChain` objects, arrays of `BubbleChain`, or arrays of `Bubble`.
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

/// Build a single JSON line (`BubbleRecord`) for a bubble by:
/// 1) building an induced context, compacting to superedges,
/// 2) splitting at endpoints,
/// 3) enumerating candidate paths with orientation-safe assembly,
/// 4) labeling via reference walk (MUE optional), then fallback adaptive flanks + KMP,
/// 5) exporting preset greedy flanks (sequence lists and node features),
/// 6) serializing **ID-free** features.
///
/// Returns `None` if bubble endpoints are malformed.
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
) -> Option<String> {
    if bubble.ends.len() != 2 { return None; }
    let start: u128 = parse_u128(&bubble.ends[0]);
    let end: u128   = parse_u128(&bubble.ends[1]);

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

    // --- Labeling cascade ---
    let mut label_opt: Option<usize> = None;
    let mut label_reason: Option<String> = None;

    if !seqs_vec.is_empty() {
        // 1) Reference-anchored walk (with MUE if endpoints are repetitive)
        let (lab1, why1) = label_with_ref_walk(
            &seqs_vec,
            &segments[&start].seq,
            &segments[&end].seq,
            reftext,
            k,
            /*try_mue=*/true
        );
        if lab1.is_some() {
            label_opt = lab1;
            label_reason = Some(why1);
        } else {
            // 2) Fallback: adaptive flanks + KMP (both strands)
            let (lab2, _bp) = label_with_adaptive_flanks(
                &seqs_vec, start, end, succ, pred, segments, k, reftext
            );
            if lab2.is_some() {
                label_opt = lab2;
                label_reason = Some("adaptive_kmp".into());
            } else {
                label_reason = Some(why1); // keep "no_ref_match"/"ambiguous_ref_multi_hit"
            }
        }
    }

    // ------- Preset flanks to export: greedily take exactly N neighbors if available -------
    let upstream_ids: Vec<u128> = collect_flanks_greedy(
        start, /*forward=*/false, succ, pred, segments, edge_cov, flank_kmers
    );
    let downstream_ids: Vec<u128> = collect_flanks_greedy(
        end, /*forward=*/true, succ, pred, segments, edge_cov, flank_kmers
    );

    // Build JSON record (ID-free)
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
        paths: paths_vec.clone(),
        label_path: label_opt,
        label_reason,
    };

    // Supernodes (as sequences)
    for &nid in &supernodes0 {
        if let Some(seg) = segments.get(&nid) {
            rec.nodes.push(NodeFeature { seq: seg.seq.clone(), cov: seg.cov });
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

/// Program entry point.
///
/// - Loads GFA graph and reference FASTA.
/// - Streams bubbles from JSON/JSONL (various shapes supported).
/// - Processes in **parallel batches** with progress bar.
/// - Writes one JSON line per bubble to `<out.jsonl>`.
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
    let radius_context: usize = 5;  // nodes to include around bubble for compaction
    let flank_kmers: usize = 5;     // preset number of upstream/downstream kmers to export
    let max_paths: usize = 256;     // candidate cap per bubble
    const BATCH: usize = 1000;      // bubbles per parallel batch (tweak for memory/CPU)

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
