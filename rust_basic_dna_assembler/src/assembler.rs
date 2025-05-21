use crate::kmer_graph::KmerGraph;
use std::collections::HashSet;

pub struct Assembler<'g> {
    pub graph: &'g KmerGraph,
}

impl<'g> Assembler<'g> {
    pub fn new(graph: &'g KmerGraph) -> Self {
        Self { graph }
    }

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
                    while self.graph.in_degree[&curr] == 1
                        && self.graph.out_degree[&curr] == 1
                    {
                        let &nxt = self.graph.edges[&curr].iter().next().unwrap();
                        visited.insert((curr, nxt));
                        curr = nxt;
                        path.push(curr);
                    }

                    let mut seq = self.graph.decode(path[0]);
                    for &code in &path[1..] {
                        let kmer = self.graph.decode(code);
                        seq.push_str(&kmer[kmer.len() - 1..]);
                    }
                    contigs.push(seq);
                }
            }
        }

        contigs
    }
}
