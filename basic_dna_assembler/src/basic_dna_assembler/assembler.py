from src.basic_dna_assembler.kmer_graph import KmerGraph
from typing import List, Tuple, Set

class Assembler:
    def __init__(self, graph: KmerGraph):
        self.graph = graph

    def assemble_contigs(self) -> List[str]:
        contigs: List[str] = []
        visited_edges: Set[Tuple[int, int]] = set()
        edges = self.graph.edges
        in_deg = self.graph.in_degree
        out_deg = self.graph.out_degree

        for u in list(edges):
            if in_deg.get(u, 0) != 1 or out_deg.get(u, 0) != 1:
                for v in edges.get(u, ()):
                    if (u, v) in visited_edges:
                        continue
                    path = [u]
                    curr = v
                    visited_edges.add((u, v))
                    while True:
                        path.append(curr)
                        if in_deg.get(curr, 0) == 1 and out_deg.get(curr, 0) == 1:
                            nbrs = edges.get(curr, ())
                            if not nbrs:
                                break
                            nxt = next(iter(nbrs))
                            visited_edges.add((curr, nxt))
                            curr = nxt
                            continue
                        break
                    seq = self.graph._decode(path[0])
                    for code in path[1:]:
                        seq += self.graph._decode(code)[-1]
                    contigs.append(seq)
        return contigs
