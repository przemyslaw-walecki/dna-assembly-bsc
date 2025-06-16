//! DNA contig assembler module.
//!
//! Provides functionality for assembling contigs from a de Bruijn graph
//! constructed from k-mers. This module defines the `Assembler` struct,
//! which identifies maximal non-branching paths in the graph and reconstructs
//! contiguous DNA sequences from them.
//!
//! # Overview
//!
//! The assembler starts contig construction at nodes with in-degree or out-degree
//! not equal to 1, and extends paths through 1-in-1-out nodes until a branch or
//! dead-end is encountered. The resulting sequences represent assembled contigs.
//!
//! # Dependencies
//!
//! - `KmerGraph`: structure representing the de Bruijn graph.
//! - `HashSet`: used to track visited edges and prevent revisiting paths.

use crate::kmer_graph::KmerGraph;
use std::collections::HashSet;


/// Struct representing a contig assembler using a de Bruijn graph.
///
/// Holds a reference to a precomputed k-mer graph and provides methods to generate contigs.
pub struct Assembler<'g> {
    /// Reference to a KmerGraph used for assembly.
    pub graph: &'g KmerGraph,
}

impl<'g> Assembler<'g> {
    /// Constructs a new assembler from a given `KmerGraph`.
    ///
    /// # Arguments
    /// * `graph` - Reference to the graph from which contigs will be assembled.
    ///
    /// # Returns
    /// A new instance of `Assembler`.
    pub fn new(graph: &'g KmerGraph) -> Self {
        Self { graph }
    }

    /// Assembles contigs from the k-mer graph by walking maximal non-branching paths.
    ///
    /// A contig is formed by starting at nodes that are not 1-in-1-out and
    /// walking through nodes with in-degree = 1 and out-degree = 1 until a branch or end is reached.
    ///
    /// # Returns
    /// A vector of contigs, each represented as a `String` DNA sequence.
    pub fn assemble_contigs(&self) -> Vec<String> {
        let mut contigs = Vec::new();
        let mut visited = HashSet::new();

        for (&u, succs) in &self.graph.edges {
            let indeg = *self.graph.in_degree.get(&u).unwrap_or(&0);
            let outdeg = *self.graph.out_degree.get(&u).unwrap_or(&0);

            // Start contig if node is a branch point or an end
            if indeg != 1 || outdeg != 1 {
                for &v in succs {
                    if !visited.insert((u, v)) {
                        continue;
                    }

                    // Walk forward while node is 1-in-1-out
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

                    // Reconstruct sequence from path of k-mers
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
