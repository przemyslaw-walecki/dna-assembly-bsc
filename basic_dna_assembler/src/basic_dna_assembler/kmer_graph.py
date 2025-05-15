from collections import defaultdict, deque
from typing import List

class KmerGraph:
    def __init__(self, k: int):
        self.k = k
        self.edges = defaultdict(list)
        self.in_degree = defaultdict(int)
        self.out_degree = defaultdict(int)

    def build(self, reads: List[str]):
        for read in reads:
            for i in range(len(read) - self.k):
                kmer1 = read[i:i+self.k]
                kmer2 = read[i+1:i+self.k+1]
                self.edges[kmer1].append(kmer2)
                self.out_degree[kmer1] += 1
                self.in_degree[kmer2] += 1

    def remove_dead_ends(self, max_depth: int = 5):
        to_remove = set()
        for node in list(self.graph.keys()):
            if len(self.graph[node]) == 0:
                path = [node]
                current = node
                for _ in range(max_depth):
                    preds = [n for n, succs in self.graph.items() if current in succs]
                    if len(preds) != 1:
                        break
                    prev = preds[0]
                    if len(self.graph[prev]) > 1:
                        break
                    path.append(prev)
                    current = prev
                if len(path) <= max_depth:
                    to_remove.update(path)

        for node in to_remove:
            self.graph.pop(node, None)
            for succs in self.graph.values():
                succs.discard(node)
    
    from collections import deque

    def remove_bubbles(self, max_depth=100):
        visited_pairs = set()
        to_remove = set()

        for source in list(self.graph.keys()):
            if len(self.graph[source]) < 2:
                continue

            paths = self._find_bubble_paths(source, max_depth)
            for path1, path2 in paths:
                if not path1 or not path2:
                    continue
                if (tuple(path1), tuple(path2)) in visited_pairs:
                    continue
                visited_pairs.add((tuple(path1), tuple(path2)))
                visited_pairs.add((tuple(path2), tuple(path1)))

                weak_path = self._choose_weaker_path(path1, path2)
                to_remove.update(weak_path)

        for node in to_remove:
            self.graph.pop(node, None)
            for succs in self.graph.values():
                succs.discard(node)

    def _find_bubble_paths(self, start, max_depth):
        paths = []
        queue = deque([(start, [start])])
        seen = {}

        while queue:
            node, path = queue.popleft()
            if len(path) > max_depth:
                continue

            for neighbor in self.graph.get(node, []):
                if neighbor in path:
                    continue
                new_path = path + [neighbor]
                if neighbor in seen:
                    paths.append((seen[neighbor], new_path))
                else:
                    seen[neighbor] = new_path
                    queue.append((neighbor, new_path))
        return paths

    def _choose_weaker_path(self, path1, path2):
        return path1 if len(path1) < len(path2) else path2

