use fxhash::FxHashMap as HashMap;
use fxhash::FxHashSet as HashSet;
use std::collections::VecDeque;

pub struct KmerGraph {
    pub k: usize,
    pub edges: HashMap<u64, HashSet<u64>>,
    pub in_degree: HashMap<u64, usize>,
    pub out_degree: HashMap<u64, usize>,
    pub edge_counts: HashMap<(u64, u64), usize>,
}

impl KmerGraph {
    /// Create a new empty graph for k-mers of size `k`.
    pub fn new(k: usize) -> Self {
        Self {
            k,
            edges: HashMap::default(),
            in_degree: HashMap::default(),
            out_degree: HashMap::default(),
            edge_counts: HashMap::default(),
        }
    }

    /// Bit-pack a k-mer string into a u64.
    pub fn encode(&self, s: &str) -> u64 {
        let mut c: u64 = 0;
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

    /// Decode a packed k-mer code back into a String.
    pub fn decode(&self, mut code: u64) -> String {
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

    /// Build the de Bruijn graph by counting k->k+1 transitions.
    pub fn build(&mut self, reads: &[String]) {
        for r in reads {
            if r.len() < self.k + 1 { continue; }
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

    /// Remove edges seen fewer than `threshold` times.
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

    /// Remove dead‐end tips up to `max_depth`.
    pub fn remove_dead_ends(&mut self, max_depth: Option<usize>) {
        let max_depth = max_depth.unwrap_or(self.k);
    
        for depth in 1..=max_depth {
            let mut pred: HashMap<u64, HashSet<u64>> = HashMap::default();
            for (&u, succs) in &self.edges {
                for &v in succs {
                    pred.entry(v).or_default().insert(u);
                }
            }
    
            let tips: Vec<u64> = self
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
                        .collect::<Vec<u64>>();
    
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
    }

    /// Remove bubbles, keeping only the longer of two distinct paths up to `max_depth`.
    pub fn remove_bubbles(&mut self, max_depth: usize) {
        let mut pred: HashMap<u64, HashSet<u64>> = HashMap::default();
        for (&u, succs) in &self.edges {
            for &v in succs {
                pred.entry(v).or_default().insert(u);
            }
        }

        let mut to_drop: HashSet<u64> = HashSet::default();
        let mut seen_pairs: HashSet<(Vec<u64>, Vec<u64>)> = HashSet::default();

        for (&src, succs) in &self.edges {
            if succs.len() < 2 {
                continue;
            }
            let mut seen: HashMap<u64, Vec<u64>> = HashMap::default();
            let mut q: VecDeque<Vec<u64>> = VecDeque::new();
            q.push_back(vec![src]);

            while let Some(path) = q.pop_front() {
                if path.len() > max_depth {
                    continue;
                }
                let n = *path.last().unwrap();
                for &nxt in self.edges.get(&n).into_iter().flat_map(|s| s.iter()) {
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
    }

    /// Total number of edges in the graph.
    pub fn edge_count(&self) -> usize {
        self.edges.values().map(|s| s.len()).sum()
    }
}
