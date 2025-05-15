class Assembler:
    def __init__(self, graph):
        self.graph = graph

    def assemble_contigs(self):
        contigs = []
        visited = set()

        def walk(node):
            contig = node
            while node in self.graph.edges and len(self.graph.edges[node]) == 1:
                next_node = self.graph.edges[node][0]
                contig += next_node[-1]
                node = next_node
            return contig

        for node in self.graph.edges:
            if self.graph.in_degree[node] != 1 or self.graph.out_degree[node] != 1:
                for next_node in self.graph.edges[node]:
                    if next_node not in visited:
                        contig = walk(next_node)
                        visited.add(next_node)
                        contigs.append(contig)
        return contigs
