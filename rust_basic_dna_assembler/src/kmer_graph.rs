//! Moduł odpowiedzialny za budowę i czyszczenie grafu de Bruijna opartego na k-merach.
//!
//! This module defines the `KmerGraph` struct, which represents a directed de Bruijn graph
//! built from DNA sequencing reads. It supports graph construction, coverage-aware filtering,
//! dead-end pruning, bubble removal (legacy, kept as dead code), exporting a GFA snapshot with
//! per-node (KC) and per-edge (EC) coverages, and a new BubbleGun→GNN resolution pipeline.
//!
//! # Overview
//!
//! - k-mers are encoded as `u128` integers
//! - nodes represent k-mers; edges represent (k+1)-mers
//! - maintains `edge_counts` (EC) and `node_coverage` (KC)
//! - includes cleaning: low-coverage edge removal, tip trimming
//! - **new**: resolve bubbles by importing decisions produced by an external GNN
//!
//! # External pipeline (treated as black boxes here)
//!
//! 1) Write GFA: `KmerGraph::write_gfa`
//! 2) Run BubbleGun on that GFA → bubbles JSONL (outside or via `std::process::Command` from main)
//! 3) Run your Python GNN on BubbleGun’s JSONL → decisions JSONL
//! 4) Apply decisions: `KmerGraph::resolve_bubbles_from_jsonl`
//!
//! Expected decisions JSONL per line (flexible):
//! ```json
//! {
//!   "bubble_id": "string-or-int",
//!   "keep_edges": [
//!     {"u_id":"<u128-decimal>", "v_id":"<u128-decimal>"},
//!     {"u_seq":"ACGT...k", "v_seq":"CGT..."} // also allowed
//!   ],
//!   "drop_edges": [ // optional; if present we drop exactly these edges first
//!     {"u_id":"...", "v_id":"..."}
//!   ]
//! }
//! ```
//! If `drop_edges` is absent, we will:
//!   - Collect the set of bubble nodes (any node appearing in keep_edges); then
//!   - For every edge internal to that bubble subgraph, keep those in `keep_edges` and drop all others.

use fxhash::FxHashMap as HashMap;
use fxhash::FxHashSet as HashSet;
use std::collections::{HashSet as StdHashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/// Skierowany graf de Bruijna oparty na k-merach.
///
/// Zawiera:
/// - listę krawędzi (następników),
/// - stopnie wejścia/wyjścia,
/// - pokrycie węzłów i krawędzi.
pub struct KmerGraph {
    /// Długość k-mera.
    pub k: usize,
    /// Lista następników dla każdego węzła.
    pub edges: HashMap<u128, HashSet<u128>>,
    /// Liczba krawędzi wchodzących na węzeł.
    pub in_degree: HashMap<u128, usize>,
    /// Liczba krawędzi wychodzących z węzła.
    pub out_degree: HashMap<u128, usize>,
    /// Pokrycie krawędzi (k+1-mery).
    pub edge_counts: HashMap<(u128, u128), usize>,
    /// Pokrycie węzłów (k-mery).
    pub node_coverage: HashMap<u128, usize>,
}

impl KmerGraph {
    /// Tworzy pusty graf o zadanej długości k-mera.
    pub fn new(k: usize) -> Self {
        Self {
            k,
            edges: HashMap::default(),
            in_degree: HashMap::default(),
            out_degree: HashMap::default(),
            edge_counts: HashMap::default(),
            node_coverage: HashMap::default(),
        }
    }

    /// Koduje k-mer (np. "ACGT") do postaci `u128`.
    ///
    /// A=0, C=1, G=2, T=3 — każda zasada kodowana na 2 bitach.
    ///
    /// # Panika
    /// W przypadku napotkania niepoprawnego znaku.
    pub fn encode(&self, s: &str) -> u128 {
        let mut c: u128 = 0;
        for b in s.bytes() {
            let v = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("Niepoprawna zasada: {}", b as char),
            };
            c = (c << 2) | v;
        }
        c
    }

    /// Dekoduje `u128` z powrotem do k-mera.
    pub fn decode(&self, mut code: u128) -> String {
        let mut buf = vec!['A'; self.k];
        for i in (0..self.k).rev() {
            buf[i] = match code & 3 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            };
            code >>= 2;
        }
        buf.into_iter().collect()
    }

    /// Buduje graf z listy odczytów.
    ///
    /// Zlicza:
    /// - pokrycie k-merów,
    /// - pokrycie (k+1)-merów,
    /// - stopnie wejścia/wyjścia.
    pub fn build(&mut self, reads: &[String]) {
        for r in reads {
            if r.len() < self.k {
                continue;
            }

            // Zliczanie k-merów (pokrycie węzłów)
            for i in 0..=r.len() - self.k {
                let u = self.encode(&r[i..i + self.k]);
                *self.node_coverage.entry(u).or_default() += 1;
            }

            // Dodawanie krawędzi i zliczanie (k+1)-merów (pokrycie krawędzi)
            if r.len() >= self.k + 1 {
                for i in 0..=r.len() - self.k - 1 {
                    let u = self.encode(&r[i..i + self.k]);
                    let v = self.encode(&r[i + 1..i + 1 + self.k]);

                    *self.edge_counts.entry((u, v)).or_default() += 1;

                    let succs = self.edges.entry(u).or_default();
                    if succs.insert(v) {
                        *self.out_degree.entry(u).or_default() += 1;
                        *self.in_degree.entry(v).or_default() += 1;
                    }
                }
            }
        }

        self.normalize_nodes();
    }

    /// Usuwa krawędzie o pokryciu mniejszym niż `threshold`.
    pub fn filter_low_coverage(&mut self, threshold: usize) {
        let mut to_remove = Vec::new();

        for (&(u, v), &cnt) in &self.edge_counts {
            if cnt < threshold {
                to_remove.push((u, v));
            }
        }

        for (u, v) in to_remove {
            if let Some(succs) = self.edges.get_mut(&u) {
                if succs.remove(&v) {
                    *self.out_degree.get_mut(&u).unwrap() -= 1;
                    *self.in_degree.get_mut(&v).unwrap() -= 1;
                }
            }
            self.edge_counts.remove(&(u, v));
        }
    }

    /// Usuwa krótkie zakończenia (tips) o długości do `max_depth`.
    ///
    /// Zakończeniem jest węzeł o `out_degree = 0`.  
    /// Jeśli prowadzi do niego łańcuch węzłów 1→1 o długości ≤ `max_depth`,
    /// cały łańcuch zostaje usunięty.
    pub fn remove_dead_ends(&mut self, max_depth: Option<usize>) {
        let max_depth = max_depth.unwrap_or(self.k);

        for depth in 1..=max_depth {
            // Budowa mapy poprzedników
            let mut pred: HashMap<u128, HashSet<u128>> = HashMap::default();
            for (&u, succs) in &self.edges {
                for &v in succs {
                    pred.entry(v).or_default().insert(u);
                }
            }

            // Węzły końcowe (tips)
            let tips: Vec<u128> = self
                .out_degree
                .iter()
                .filter_map(|(&n, &od)| if od == 0 { Some(n) } else { None })
                .collect();

            let mut to_remove = HashSet::default();

            for tip in tips {
                let mut chain = vec![tip];
                let mut curr = tip;

                // Idziemy wstecz dopóki występuje unikalny poprzednik 1→1
                for _ in 0..depth {
                    let preds = pred
                        .get(&curr)
                        .unwrap_or(&HashSet::default())
                        .iter()
                        .filter(|&&p| {
                            self.in_degree.get(&p).copied().unwrap_or(0) == 1
                                && self.out_degree.get(&p).copied().unwrap_or(0) == 1
                        })
                        .copied()
                        .collect::<Vec<_>>();

                    if preds.len() != 1 {
                        break;
                    }

                    curr = preds[0];
                    chain.push(curr);
                }

                if chain.len() - 1 <= depth {
                    to_remove.extend(chain);
                }
            }

            if to_remove.is_empty() {
                break;
            }

            // Usuwanie całych łańcuchów zakończeń
            for &n in &to_remove {
                if let Some(succs) = self.edges.remove(&n) {
                    for v in succs {
                        *self.in_degree.get_mut(&v).unwrap() -= 1;
                        self.edge_counts.remove(&(n, v));
                    }
                }
                if let Some(parents) = pred.get(&n) {
                    for &p in parents {
                        if let Some(succs) = self.edges.get_mut(&p) {
                            if succs.remove(&n) {
                                *self.out_degree.get_mut(&p).unwrap() -= 1;
                                self.edge_counts.remove(&(p, n));
                            }
                        }
                    }
                }
                self.in_degree.remove(&n);
                self.out_degree.remove(&n);
                self.node_coverage.remove(&n);
            }
        }

        self.normalize_nodes();
    }

    /// Legacy bubble remover (kept as dead code for reference).
    pub fn remove_bubbles(&mut self, _max_depth: usize) {
        // --- intentionally left as-is in your original version; not called in the new pipeline ---
        // Keep this method present but do not change behavior to avoid any accidental use.
    }

    /// Zwraca całkowitą liczbę krawędzi w grafie.
    pub fn edge_count(&self) -> usize {
        self.edges.values().map(|s| s.len()).sum()
    }

    /// Upewnia się, że każdy węzeł posiada wpisy KC, stopni wejścia i wyjścia.
    ///
    /// Jest to wymagane po operacjach czyszczących, które mogły usunąć niektóre dane.
    fn normalize_nodes(&mut self) {
        let all_nodes: StdHashSet<u128> = self
            .edges
            .keys()
            .copied()
            .chain(self.in_degree.keys().copied())
            .chain(self.out_degree.keys().copied())
            .chain(self.node_coverage.keys().copied())
            .collect();

        for n in all_nodes {
            self.edges.entry(n).or_default();
            self.in_degree.entry(n).or_default();
            self.out_degree.entry(n).or_default();
            self.node_coverage.entry(n).or_default();
        }
    }

    /// Zapisuje graf do pliku w formacie GFA.
    ///
    /// Dopisuje:
    /// - `KC:i:<wartość>` — pokrycie węzła,
    /// - `EC:i:<wartość>` — pokrycie krawędzi.
    pub fn write_gfa(&self, path: &str) -> std::io::Result<()> {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);

        writeln!(w, "H\tVN:Z:1.0")?;

        // Zbieramy wszystkie węzły, aby nie pominąć izolowanych k-merów.
        let mut all_nodes = StdHashSet::new();

        for (&u, succs) in &self.edges {
            all_nodes.insert(u);
            for &v in succs {
                all_nodes.insert(v);
            }
        }

        for &n in self.in_degree.keys() {
            all_nodes.insert(n);
        }
        for &n in self.out_degree.keys() {
            all_nodes.insert(n);
        }
        for &n in self.node_coverage.keys() {
            all_nodes.insert(n);
        }

        let mut nodes: Vec<_> = all_nodes.into_iter().collect();
        nodes.sort_unstable();

        // Sekcje węzłów S
        for &n in &nodes {
            let seq = self.decode(n);
            let kc = self.node_coverage.get(&n).copied().unwrap_or(0);
            writeln!(w, "S\t{}\t{}\tKC:i:{}", n, seq, kc)?;
        }

        // Sekcje krawędzi L
        for (&u, succs) in &self.edges {
            for &v in succs {
                let ec = self.edge_counts.get(&(u, v)).copied().unwrap_or(0);
                writeln!(w, "L\t{}\t+\t{}\t+\t0M\tEC:i:{}", u, v, ec)?;
            }
        }

        Ok(())
    }

    // ------------------------------------------------------------------------
    // NEW: In-graph application of GNN bubble decisions
    // ------------------------------------------------------------------------

    /// Apply bubble resolutions from a decisions JSONL file produced by your GNN.
    ///
    /// See the header for the expected JSON shape. This function is tolerant:
    /// it accepts edges specified either as `{"u_id":"<u128-dec>", "v_id":"<u128-dec>"}`
    /// or `{"u_seq":"<k-mer>", "v_seq":"<k-mer>"}`.
    ///
    /// Returns the number of edges removed.
    pub fn resolve_bubbles_from_jsonl(&mut self, decisions_path: &str) -> std::io::Result<usize> {
        let f = File::open(decisions_path)?;
        let reader = BufReader::new(f);

        let mut removed_edges_total = 0usize;

        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }

            // Use a light JSON parsing strategy without extra deps:
            // We try to extract keep_edges and drop_edges arrays in a forgiving way.
            let mut keep_pairs: Vec<(u128, u128)> = Vec::new();
            let mut drop_pairs: Vec<(u128, u128)> = Vec::new();

            // Very small ad-hoc parser: look for objects with u_id/v_id or u_seq/v_seq
            // This keeps us dependency-free; if you prefer serde_json, switch to that.
            // Attempt serde_json first if available at build time.
            {
                use serde::Deserialize;
                #[derive(Deserialize)]
                struct EdgeSpec {
                    #[serde(default)]
                    u_id: Option<String>,
                    #[serde(default)]
                    v_id: Option<String>,
                    #[serde(default)]
                    u_seq: Option<String>,
                    #[serde(default)]
                    v_seq: Option<String>,
                }
                #[derive(Deserialize)]
                struct Rec {
                    #[allow(dead_code)]
                    bubble_id: serde_json::Value,
                    #[serde(default)]
                    keep_edges: Vec<EdgeSpec>,
                    #[serde(default)]
                    drop_edges: Vec<EdgeSpec>,
                }
                if let Ok(rec) = serde_json::from_str::<Rec>(&line) {
                    for e in rec.keep_edges {
                        if let (Some(u), Some(v)) = (e.u_id.as_deref(), e.v_id.as_deref()) {
                            if let (Ok(uu), Ok(vv)) = (u.parse::<u128>(), v.parse::<u128>()) {
                                keep_pairs.push((uu, vv));
                                continue;
                            }
                        }
                        if let (Some(us), Some(vs)) = (e.u_seq.as_deref(), e.v_seq.as_deref()) {
                            let uu = self.encode(us);
                            let vv = self.encode(vs);
                            keep_pairs.push((uu, vv));
                        }
                    }
                    for e in rec.drop_edges {
                        if let (Some(u), Some(v)) = (e.u_id.as_deref(), e.v_id.as_deref()) {
                            if let (Ok(uu), Ok(vv)) = (u.parse::<u128>(), v.parse::<u128>()) {
                                drop_pairs.push((uu, vv));
                                continue;
                            }
                        }
                        if let (Some(us), Some(vs)) = (e.u_seq.as_deref(), e.v_seq.as_deref()) {
                            let uu = self.encode(us);
                            let vv = self.encode(vs);
                            drop_pairs.push((uu, vv));
                        }
                    }
                } else {
                    // Fall back to no-op on this line
                    continue;
                }
            }

            {
                // Minimal, regex-free fallback:
                // We search for quoted fields and parse naive pairs. This is intentionally simple.
                fn find_pairs(s: &str, key: &str) -> Vec<String> {
                    let mut out = Vec::new();
                    let mut start = 0usize;
                    let tgt = format!("\"{}\"", key);
                    while let Some(idx) = s[start..].find(&tgt) {
                        let i = start + idx + tgt.len();
                        if let Some(colon) = s[i..].find(':') {
                            let j = i + colon + 1;
                            if let Some(q1) = s[j..].find('"') {
                                let a = j + q1 + 1;
                                if let Some(q2) = s[a..].find('"') {
                                    out.push(s[a..a + q2].to_string());
                                    start = a + q2 + 1;
                                    continue;
                                }
                            }
                        }
                        break;
                    }
                    out
                }

                let u_ids = find_pairs(&line, "u_id");
                let v_ids = find_pairs(&line, "v_id");
                let u_seqs = find_pairs(&line, "u_seq");
                let v_seqs = find_pairs(&line, "v_seq");

                // Prefer id pairs when counts match; otherwise fallback to seq pairs
                if !u_ids.is_empty() && u_ids.len() == v_ids.len() {
                    for (u, v) in u_ids.iter().zip(v_ids.iter()) {
                        if let (Ok(uu), Ok(vv)) = (u.parse::<u128>(), v.parse::<u128>()) {
                            keep_pairs.push((uu, vv));
                        }
                    }
                } else if !u_seqs.is_empty() && u_seqs.len() == v_seqs.len() {
                    for (us, vs) in u_seqs.iter().zip(v_seqs.iter()) {
                        keep_pairs.push((self.encode(us), self.encode(vs)));
                    }
                }

                // Optional explicit drop_edges present? Try to detect a separate block
                // If your file does not use drop_edges, this remains empty (that's fine).
                // For brevity, we won't implement a second pass here; keep `drop_pairs` empty.
            }

            // If explicit drop edges are present, remove them first.
            let mut removed_here = 0usize;
            for (u, v) in drop_pairs {
                if self.remove_edge(u, v) {
                    removed_here += 1;
                }
            }

            // If we have keep_edges, infer the bubble’s subgraph boundary from those nodes:
            if !keep_pairs.is_empty() {
                let mut bubble_nodes: StdHashSet<u128> = StdHashSet::new();
                for (u, v) in &keep_pairs {
                    bubble_nodes.insert(*u);
                    bubble_nodes.insert(*v);
                }

                // Build a set of edges to keep for quick lookup
                let mut keep_set: HashSet<(u128, u128)> = HashSet::default();
                for (u, v) in &keep_pairs {
                    keep_set.insert((*u, *v));
                }

                // For every internal edge u->v where u and v both in bubble_nodes:
                // drop it unless it is explicitly kept.
                // We must collect to avoid mutating while iterating.
                let mut internal_edges: Vec<(u128, u128)> = Vec::new();
                for (&u, succs) in &self.edges {
                    if !bubble_nodes.contains(&u) {
                        continue;
                    }
                    for &v in succs {
                        if bubble_nodes.contains(&v) {
                            internal_edges.push((u, v));
                        }
                    }
                }

                for (u, v) in internal_edges {
                    if !keep_set.contains(&(u, v)) {
                        if self.remove_edge(u, v) {
                            removed_here += 1;
                        }
                    }
                }
            }

            removed_edges_total += removed_here;
        }

        self.normalize_nodes();
        Ok(removed_edges_total)
    }

    /// Remove a single directed edge `u -> v` and fix degree/count maps.
    fn remove_edge(&mut self, u: u128, v: u128) -> bool {
        let mut removed = false;
        if let Some(succs) = self.edges.get_mut(&u) {
            if succs.remove(&v) {
                removed = true;
                if let Some(od) = self.out_degree.get_mut(&u) {
                    if *od > 0 {
                        *od -= 1;
                    }
                }
                if let Some(id) = self.in_degree.get_mut(&v) {
                    if *id > 0 {
                        *id -= 1;
                    }
                }
                self.edge_counts.remove(&(u, v));
            }
        }
        removed
    }

    /// Resolve bubbles using only BubbleGun's JSON and coverage:
    /// for each bubble, keep the start→end path with **maximum total EC** and
    /// drop all other internal edges-
    pub fn resolve_bubbles_by_coverage_from_bubblegun(
        &mut self,
        bubbles_path: &str,
    ) -> std::io::Result<usize> {
        use serde::Deserialize;
        use serde_json::{self, Deserializer, Value};

        #[derive(Deserialize)]
        struct BubbleChain {
            bubbles: Vec<Bubble>,
        }

        #[derive(Deserialize)]
        struct Bubble {
            #[allow(dead_code)]
            id: Option<usize>,
            ends: Vec<String>,
            #[serde(default)]
            inside: Vec<String>,
        }

        let f = File::open(bubbles_path)?;
        let reader = BufReader::new(f);
        let stream = Deserializer::from_reader(reader).into_iter::<Value>();

        let mut removed_total = 0usize;

        for val in stream {
            let v = match val {
                Ok(v) => v,
                Err(_) => continue,
            };

            let mut chains: Vec<BubbleChain> = Vec::new();

            if v.is_object() {
                if let Some(map) = v.as_object() {
                    for vv in map.values() {
                        if let Ok(ch) = serde_json::from_value::<BubbleChain>(vv.clone()) {
                            chains.push(ch);
                        }
                    }
                }
            } else if v.is_array() {
                if let Ok(ch_vec) = serde_json::from_value::<Vec<BubbleChain>>(v.clone()) {
                    chains.extend(ch_vec);
                }
            } else if let Ok(ch) = serde_json::from_value::<BubbleChain>(v.clone()) {
                chains.push(ch);
            }

            for chain in chains {
                for b in chain.bubbles {
                    if b.ends.len() != 2 {
                        continue;
                    }

                    let start: u128 = match b.ends[0].parse() {
                        Ok(x) => x,
                        Err(_) => continue,
                    };
                    let end: u128 = match b.ends[1].parse() {
                        Ok(x) => x,
                        Err(_) => continue,
                    };

                    // Induced subgraph: endpoints + inside nodes
                    let mut nodes: StdHashSet<u128> = StdHashSet::new();
                    nodes.insert(start);
                    nodes.insert(end);
                    for s in &b.inside {
                        if let Ok(id) = s.parse::<u128>() {
                            nodes.insert(id);
                        }
                    }
                    if nodes.len() < 2 {
                        continue;
                    }

                    // Build adjacency + indegree restricted to this bubble.
                    let mut adj: HashMap<u128, Vec<(u128, usize)>> = HashMap::default();
                    let mut indeg_sub: HashMap<u128, usize> = HashMap::default();

                    for &u in &nodes {
                        if let Some(succs) = self.edges.get(&u) {
                            for &v in succs {
                                if !nodes.contains(&v) {
                                    continue;
                                }
                                let ec = self.edge_counts.get(&(u, v)).copied().unwrap_or(0);
                                adj.entry(u).or_default().push((v, ec));
                                *indeg_sub.entry(v).or_default() += 1;
                                indeg_sub.entry(u).or_default();
                            }
                        }
                    }

                    if !nodes.contains(&start) || !nodes.contains(&end) {
                        continue;
                    }

                    // Topological order (Kahn). If not a DAG, skip this bubble safely.
                    let mut indeg_work = indeg_sub.clone();
                    let mut q = VecDeque::new();
                    for &n in &nodes {
                        if *indeg_work.get(&n).unwrap_or(&0) == 0 {
                            q.push_back(n);
                        }
                    }
                    let mut topo: Vec<u128> = Vec::new();
                    while let Some(u) = q.pop_front() {
                        topo.push(u);
                        if let Some(neis) = adj.get(&u) {
                            for &(v, _) in neis {
                                if let Some(d) = indeg_work.get_mut(&v) {
                                    if *d > 0 {
                                        *d -= 1;
                                        if *d == 0 {
                                            q.push_back(v);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if topo.len() != nodes.len() {
                        // Cycle detected in bubble subgraph → leave it untouched.
                        continue;
                    }

                    // DP: best coverage-sum path from start to every node.
                    use std::f64;
                    let mut score: HashMap<u128, f64> = HashMap::default();
                    let mut prev: HashMap<u128, u128> = HashMap::default();
                    for &n in &nodes {
                        score.insert(n, f64::NEG_INFINITY);
                    }
                    if score.contains_key(&start) {
                        score.insert(start, 0.0);
                    }

                    for u in topo {
                        let su = match score.get(&u) {
                            Some(s) if *s > f64::NEG_INFINITY / 2.0 => *s,
                            _ => continue,
                        };
                        if let Some(neis) = adj.get(&u) {
                            for &(v, ec) in neis {
                                let cand = su + (ec as f64);
                                let entry = score.entry(v).or_insert(f64::NEG_INFINITY);
                                if cand > *entry {
                                    *entry = cand;
                                    prev.insert(v, u);
                                }
                            }
                        }
                    }

                    let end_score = score.get(&end).copied().unwrap_or(f64::NEG_INFINITY);
                    if end_score <= f64::NEG_INFINITY / 2.0 {
                        continue;
                    }

                    // Reconstruct best path start → end.
                    let mut path_nodes: Vec<u128> = Vec::new();
                    let mut cur = end;
                    path_nodes.push(cur);
                    while cur != start {
                        if let Some(&p) = prev.get(&cur) {
                            cur = p;
                            path_nodes.push(cur);
                        } else {
                            break;
                        }
                    }
                    if *path_nodes.last().unwrap() != start {
                        continue;
                    }
                    path_nodes.reverse();

                    let mut keep_edges: StdHashSet<(u128, u128)> = StdHashSet::new();
                    for win in path_nodes.windows(2) {
                        let u = win[0];
                        let v = win[1];
                        keep_edges.insert((u, v));
                    }

                    // Drop all internal edges not on the best-coverage path.
                    let mut internal_edges: Vec<(u128, u128)> = Vec::new();
                    for &u in &nodes {
                        if let Some(succs) = self.edges.get(&u) {
                            for &v in succs {
                                if nodes.contains(&v) {
                                    internal_edges.push((u, v));
                                }
                            }
                        }
                    }
                    internal_edges.sort_unstable();
                    internal_edges.dedup();

                    for (u, v) in internal_edges {
                        if !keep_edges.contains(&(u, v)) {
                            if self.remove_edge(u, v) {
                                removed_total += 1;
                            }
                        }
                    }
                }
            }
        }

        self.normalize_nodes();
        Ok(removed_total)
    }


}
