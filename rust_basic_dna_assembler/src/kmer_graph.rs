use fxhash::FxHashMap as HashMap;
use fxhash::FxHashSet as HashSet;
use std::collections::{HashSet as StdHashSet, VecDeque};
use std::fs::File;
use std::io::{BufWriter, Write};

/// Directed de Bruijn graph for k-mers.
pub struct KmerGraph {
    pub k: usize,
    pub edges: HashMap<u128, HashSet<u128>>,
    pub in_degree: HashMap<u128, usize>,
    pub out_degree: HashMap<u128, usize>,
    /// (k+1)-mer counts (edge coverage)
    pub edge_counts: HashMap<(u128, u128), usize>,
    /// k-mer counts (node coverage)
    pub node_counts: HashMap<u128, usize>,
}

impl KmerGraph {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            edges: HashMap::default(),
            in_degree: HashMap::default(),
            out_degree: HashMap::default(),
            edge_counts: HashMap::default(),
            node_counts: HashMap::default(),
        }
    }

    pub fn encode(&self, s: &str) -> u128 {
        let mut c: u128 = 0;
        for b in s.bytes() {
            let v = match b {
                b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3,
                _ => panic!("Invalid base {}", b as char),
            };
            c = (c << 2) | v;
        }
        c
    }

    pub fn decode(&self, mut code: u128) -> String {
        let mut buf = vec!['A'; self.k];
        for i in (0..self.k).rev() {
            buf[i] = match code & 3 { 0 => 'A', 1 => 'C', 2 => 'G', _ => 'T' };
            code >>= 2;
        }
        buf.into_iter().collect()
    }

    /// Build graph and count both k-mers (KC) and (k+1)-mers (EC)
    pub fn build(&mut self, reads: &[String]) {
        for r in reads {
            if r.len() < self.k { continue; }

            // count k-mers (node coverage)
            for i in 0..=r.len() - self.k {
                let u = self.encode(&r[i..i + self.k]);
                *self.node_counts.entry(u).or_default() += 1;
            }

            // add edges and count (k+1)-mers (edge coverage)
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
    }

    pub fn filter_low_coverage(&mut self, threshold: usize) {
        let mut to_remove = Vec::new();
        for (&(u, v), &cnt) in &self.edge_counts {
            if cnt < threshold { to_remove.push((u, v)); }
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

    // ... (dead-end trimming / bubbles unchanged) ...

    pub fn edge_count(&self) -> usize {
        self.edges.values().map(|s| s.len()).sum()
    }

    /// Ensure degree maps include every node we know about
    fn normalize_nodes(&mut self) {
        let all_nodes: StdHashSet<u128> = self.edges.keys()
            .chain(self.in_degree.keys())
            .chain(self.out_degree.keys())
            .chain(self.node_counts.keys())   // <- include k-mer count nodes
            .copied()
            .collect();

        for n in all_nodes {
            self.edges.entry(n).or_default();
            self.in_degree.entry(n).or_default();
            self.out_degree.entry(n).or_default();
            self.node_counts.entry(n).or_default();
        }
    }

    /// Write GFA with KC (nodes) and EC (edges)
    pub fn write_gfa(&self, path: &str) -> std::io::Result<()> {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);

        writeln!(w, "H\tVN:Z:1.0")?;

        // union of everything so we don't drop isolated k-mers
        let mut all_nodes = StdHashSet::new();
        for (&u, succs) in &self.edges {
            all_nodes.insert(u);
            for &v in succs { all_nodes.insert(v); }
        }
        for &n in self.in_degree.keys() { all_nodes.insert(n); }
        for &n in self.out_degree.keys() { all_nodes.insert(n); }
        for &n in self.node_counts.keys() { all_nodes.insert(n); }

        let mut nodes: Vec<_> = all_nodes.into_iter().collect();
        nodes.sort_unstable();

        // S-lines with KC:i
        for &n in &nodes {
            let seq = self.decode(n);
            let kc = self.node_counts.get(&n).copied().unwrap_or(0);
            writeln!(w, "S\t{}\t{}\tKC:i:{}", n, seq, kc)?;
        }

        // L-lines with EC:i
        for (&u, succs) in &self.edges {
            for &v in succs {
                let ec = self.edge_counts.get(&(u, v)).copied().unwrap_or(0);
                writeln!(w, "L\t{}\t+\t{}\t+\t0M\tEC:i:{}", u, v, ec)?;
            }
        }

        Ok(())
    }
}
