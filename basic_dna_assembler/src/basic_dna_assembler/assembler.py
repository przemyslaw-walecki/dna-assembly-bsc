class Assembler:
    def __init__(self, graph):
        """
        Walk non-1->1 regions to extract contigs.
        """
        self.graph = graph

    def assemble_contigs(self):
        contigs = []
        visited = set()
        for u, succs in self.graph.edges.items():
            if self.graph.in_degree[u] != 1 or self.graph.out_degree[u] != 1:
                for v in succs:
                    if (u, v) in visited:
                        continue
                    path = [u, v]
                    visited.add((u, v))
                    curr = v
                    while (
                        self.graph.in_degree[curr] == 1
                        and self.graph.out_degree[curr] == 1
                    ):
                        nxt = next(iter(self.graph.edges[curr]))
                        visited.add((curr, nxt))
                        curr = nxt
                        path.append(curr)

                    seq = self.graph._decode(path[0])
                    for code in path[1:]:
                        seq += self.graph._decode(code)[-1]
                    contigs.append(seq)
        return contigs
