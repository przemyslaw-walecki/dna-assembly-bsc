//! K-mer de Bruijn graph construction and simplification module.
//!
//! This module defines the `KmerGraph` struct, which represents a directed de Bruijn graph
//! built from DNA sequencing reads. It includes methods for graph construction, edge filtering,
//! dead-end pruning, bubble removal, and decoding encoded k-mers.
//!
//! # Overview
//!
//! - k-mers are encoded as `u128` integers
//! - nodes represent k-mers; edges represent (k+1)-mers
//! - supports graph cleaning operations like low-coverage removal, dead-end trimming, and bubble removal
//!
//! # Dependencies
//!
//! - Uses `fxhash` for fast hashing
//! - `VecDeque` for breadth-first search during bubble detection

use fxhash::FxHashMap as HashMap;
use fxhash::FxHashSet as HashSet;
use std::collections::{HashSet as StdHashSet, VecDeque};
use std::fs::File;
use std::io::{BufWriter, Write};

/// Directed de Bruijn graph for k-mers.
pub struct KmerGraph {
    /// Length of k-mers used in the graph.
    pub k: usize,
    /// Adjacency list: maps each node to a set of successor nodes.
    pub edges: HashMap<u128, HashSet<u128>>,
    /// Maps node to in-degree count.
    pub in_degree: HashMap<u128, usize>,
    /// Maps node to out-degree count.
    pub out_degree: HashMap<u128, usize>,
    /// Number of times each edge (u, v) appeared in reads.
    pub edge_counts: HashMap<(u128, u128), usize>,
}

impl KmerGraph {
    /// Creates a new, empty `KmerGraph` with the given k-mer length.
    pub fn new(k: usize) -> Self {
        Self {
            k,
            edges: HashMap::default(),
            in_degree: HashMap::default(),
            out_degree: HashMap::default(),
            edge_counts: HashMap::default(),
        }
    }

    /// Encodes a DNA k-mer string (e.g., `"ACG"`) into a compact `u128` integer.
    ///
    /// # Panics
    /// If the string contains invalid bases.
    pub fn encode(&self, s: &str) -> u128 {
        let mut c: u128 = 0;
        for b in s.bytes() {
            let v = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("Invalid base {}", b as char),
            };
            c = (c << 2) | v;
        }
        c
    }

    /// Decodes a `u128` integer into a DNA k-mer string (e.g., `"ACG"`).
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

    /// Builds the de Bruijn graph from a list of sequencing reads.
    ///
    /// Adds edges and updates in/out degrees and edge counts.
    pub fn build(&mut self, reads: &[String]) {
        for r in reads {
            if r.len() < self.k + 1 {
                continue;
            }
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

    /// Removes edges with coverage below the given threshold.
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
        }
    }

    /// Removes short dead-end branches (tips) of the graph up to a specified depth.
    ///
    /// # Arguments
    /// * `max_depth` - maximum chain length to remove (defaults to `k` if `None`)
    pub fn remove_dead_ends(&mut self, max_depth: Option<usize>) {
        let max_depth = max_depth.unwrap_or(self.k);

        for depth in 1..=max_depth {
            let mut pred: HashMap<u128, HashSet<u128>> = HashMap::default();
            for (&u, succs) in &self.edges {
                for &v in succs {
                    pred.entry(v).or_default().insert(u);
                }
            }

            let tips: Vec<u128> = self
                .out_degree
                .iter()
                .filter_map(|(&n, &od)| if od == 0 { Some(n) } else { None })
                .collect();

            let mut to_remove = HashSet::default();
            for tip in tips {
                let mut chain = vec![tip];
                let mut curr = tip;
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
                        .collect::<Vec<u128>>();

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

            for &n in &to_remove {
                if let Some(succs) = self.edges.remove(&n) {
                    for v in succs {
                        *self.in_degree.get_mut(&v).unwrap() -= 1;
                    }
                }
                if let Some(parents) = pred.get(&n) {
                    for &p in parents {
                        if let Some(succs) = self.edges.get_mut(&p) {
                            if succs.remove(&n) {
                                *self.out_degree.get_mut(&p).unwrap() -= 1;
                            }
                        }
                    }
                }
                self.in_degree.remove(&n);
                self.out_degree.remove(&n);
            }
        }

        self.normalize_nodes();
    }

    /// Removes bubbles from the graph (redundant alternative paths).
    ///
    /// # Arguments
    /// * `max_depth` - maximum length of bubble paths to detect and collapse.
    pub fn remove_bubbles(&mut self, max_depth: usize) {
        let mut pred: HashMap<u128, HashSet<u128>> = HashMap::default();
        for (&u, succs) in &self.edges {
            for &v in succs {
                pred.entry(v).or_default().insert(u);
            }
        }

        let mut to_drop: HashSet<u128> = HashSet::default();
        let mut seen_pairs: HashSet<(Vec<u128>, Vec<u128>)> = HashSet::default();

        for (&src, succs) in &self.edges {
            if succs.len() < 2 {
                continue;
            }
            let mut seen: HashMap<u128, Vec<u128>> = HashMap::default();
            let mut q: VecDeque<Vec<u128>> = VecDeque::new();
            q.push_back(vec![src]);

            while let Some(path) = q.pop_front() {
                if path.len() > max_depth {
                    continue;
                }
                let n = *path.last().unwrap();
                if let Some(neighbors) = self.edges.get(&n) {
                    for &nxt in neighbors {
                        if path.contains(&nxt) {
                            continue;
                        }
                        let mut newp = path.clone();
                        newp.push(nxt);
                        if let Some(p2) = seen.get(&nxt) {
                            let key = (newp.clone(), p2.clone());
                            if !seen_pairs.contains(&key) {
                                seen_pairs.insert(key.clone());
                                seen_pairs.insert((p2.clone(), newp.clone()));
                                let drop_path = if newp.len() < p2.len() { &newp } else { p2 };
                                to_drop.extend(drop_path.iter().copied());
                            }
                        } else {
                            seen.insert(nxt, newp.clone());
                            q.push_back(newp);
                        }
                    }
                }
            }
        }

        for &n in &to_drop {
            if let Some(succs) = self.edges.remove(&n) {
                for v in succs {
                    *self.in_degree.get_mut(&v).unwrap() -= 1;
                }
            }
            if let Some(parents) = pred.get(&n) {
                for &p in parents {
                    if let Some(succs) = self.edges.get_mut(&p) {
                        if succs.remove(&n) {
                            *self.out_degree.get_mut(&p).unwrap() -= 1;
                        }
                    }
                }
            }
            self.in_degree.remove(&n);
            self.out_degree.remove(&n);
        }

        self.normalize_nodes();
    }

    /// Returns the number of edges in the graph.
    pub fn edge_count(&self) -> usize {
        self.edges.values().map(|s| s.len()).sum()
    }

    /// Ensures that all nodes in the graph have in/out-degree entries.
    fn normalize_nodes(&mut self) {
        let all_nodes: StdHashSet<u128> = self
            .edges
            .keys()
            .chain(self.in_degree.keys())
            .chain(self.out_degree.keys())
            .copied()
            .collect();
        for n in all_nodes {
            self.edges.entry(n).or_default();
            self.in_degree.entry(n).or_default();
            self.out_degree.entry(n).or_default();
        }
    }
    
    pub fn write_gfa(&self, path: &str) -> std::io::Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
    
        writeln!(writer, "H\tVN:Z:1.0")?;
    
        // Collect all nodes: from edges, in_degree, out_degree
        let mut all_nodes = StdHashSet::new();
        for &u in self.edges.keys() {
            all_nodes.insert(u);
            if let Some(vs) = self.edges.get(&u) {
                for &v in vs {
                    all_nodes.insert(v);
                }
            }
        }
    
        // Write segments (S lines)
        let mut all_nodes_vec: Vec<_> = all_nodes.into_iter().collect();
        all_nodes_vec.sort_unstable();
        for node in all_nodes_vec.iter() {
            let seq = self.decode(*node);
            writeln!(writer, "S\t{}\t{}", node, seq)?;
        }
    
        // Write links (L lines)
        for (&u, succs) in &self.edges {
            for &v in succs {
                writeln!(writer, "L\t{}\t+\t{}\t+\t0M", u, v)?;
            }
        }
    
        Ok(())
    }
}
