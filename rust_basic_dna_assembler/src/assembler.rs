//! DNA contig assembler module.
//!
//! Exports maximal unitigs (including pure 1-in-1-out cycles). The previous
//! assemble_contigs() is preserved as dead code and not called.

use crate::kmer_graph::KmerGraph;
use std::collections::HashSet;

pub struct Assembler<'g> {
    pub graph: &'g KmerGraph,
}

impl<'g> Assembler<'g> {
    pub fn new(graph: &'g KmerGraph) -> Self {
        Self { graph }
    }

    // --- OLD METHOD (dead code; not called) ---------------------------------
    pub fn assemble_contigs(&self) -> Vec<String> {
        let mut contigs = Vec::new();
        let mut visited = HashSet::new();

        for (&u, succs) in &self.graph.edges {
            let indeg = *self.graph.in_degree.get(&u).unwrap_or(&0);
            let outdeg = *self.graph.out_degree.get(&u).unwrap_or(&0);

            if indeg != 1 || outdeg != 1 {
                for &v in succs {
                    if !visited.insert((u, v)) {
                        continue;
                    }

                    let mut path = vec![u, v];
                    let mut curr = v;

                    // SAFE access
                    while self.deg_in(curr) == 1 && self.deg_out(curr) == 1 {
                        if let Some(nxt) = self.first_succ(curr) {
                            if !visited.insert((curr, nxt)) {
                                break;
                            }
                            curr = nxt;
                            path.push(curr);
                        } else {
                            break;
                        }
                    }

                    contigs.push(self.reconstruct(&path, false));
                }
            }
        }
        contigs
    }

    // --- NEW METHOD (unitigs incl. cycles) ----------------------------------
    pub fn assemble_unitigs(&self) -> Vec<String> {
        let mut unitigs = Vec::new();
        let mut visited: HashSet<(u128, u128)> = HashSet::new();

        // Step 1: start at branch/end nodes
        for (&u, succs) in &self.graph.edges {
            let indeg = self.deg_in(u);
            let outdeg = self.deg_out(u);

            if indeg != 1 || outdeg != 1 {
                for &v in succs {
                    if !visited.insert((u, v)) {
                        continue;
                    }

                    let mut path = vec![u, v];
                    let mut curr = v;

                    while self.deg_in(curr) == 1 && self.deg_out(curr) == 1 {
                        if let Some(nxt) = self.first_succ(curr) {
                            if !visited.insert((curr, nxt)) {
                                break;
                            }
                            curr = nxt;
                            path.push(curr);
                        } else {
                            // out_degree said 1 but adjacency missing; stop gracefully
                            break;
                        }
                    }

                    unitigs.push(self.reconstruct(&path, false));
                }
            }
        }

        // Step 2: handle pure 1-in-1-out cycles
        for (&u, succs) in &self.graph.edges {
            if self.deg_in(u) != 1 || self.deg_out(u) != 1 {
                continue;
            }

            for &v in succs {
                if visited.contains(&(u, v)) {
                    continue;
                }

                // Walk cycle starting from edge (u, v)
                let mut path = vec![u, v];
                visited.insert((u, v));
                let mut curr = v;

                loop {
                    // Leave if not 1-in-1-out (should be rare after step 1 handled branches)
                    if self.deg_in(curr) != 1 || self.deg_out(curr) != 1 {
                        // Extend linearly until a branch/end; emit as linear path
                        while self.deg_in(curr) == 1 && self.deg_out(curr) == 1 {
                            if let Some(nxt) = self.first_succ(curr) {
                                if !visited.insert((curr, nxt)) {
                                    break;
                                }
                                curr = nxt;
                                path.push(curr);
                            } else {
                                break;
                            }
                        }
                        unitigs.push(self.reconstruct(&path, false));
                        break;
                    }

                    let Some(nxt) = self.first_succ(curr) else {
                        // Inconsistent maps; emit what we have
                        unitigs.push(self.reconstruct(&path, false));
                        break;
                    };

                    // If we are about to traverse the starting edge again, the cycle is closed.
                    if curr == u && nxt == v {
                        // Push start once to mark closure for reconstruction logic
                        path.push(u);
                        unitigs.push(self.reconstruct(&path, true));
                        break;
                    }

                    if !visited.insert((curr, nxt)) {
                        // Another walk already covered the remainder; nothing to emit
                        break;
                    }

                    curr = nxt;
                    path.push(curr);
                }
            }
        }

        unitigs
    }

    // --- Helpers -------------------------------------------------------------

    #[inline]
    fn deg_in(&self, u: u128) -> usize {
        *self.graph.in_degree.get(&u).unwrap_or(&0)
    }

    #[inline]
    fn deg_out(&self, u: u128) -> usize {
        *self.graph.out_degree.get(&u).unwrap_or(&0)
    }

    #[inline]
    fn first_succ(&self, u: u128) -> Option<u128> {
        self.graph
            .edges
            .get(&u)
            .and_then(|s| s.iter().next().copied())
    }

    /// Reconstruct sequence from a path of node codes.
    /// For cycles, `path` ends with the start node duplicated once; we avoid
    /// duplicating the last base from the duplicated node.
    fn reconstruct(&self, path: &[u128], is_cycle: bool) -> String {
        if path.is_empty() {
            return String::new();
        }
        let mut seq = self.graph.decode(path[0]);
        let end = if is_cycle {
            // Last entry equals the first; skip it when appending bases
            path.len() - 1
        } else {
            path.len()
        };
        for &code in &path[1..end] {
            let kmer = self.graph.decode(code);
            // append last base
            seq.push_str(&kmer[kmer.len() - 1..]);
        }
        seq
    }
}
