from collections import defaultdict, deque
from typing import List, Set, Tuple

class KmerGraph:
    def __init__(self, k: int):
        self.k = k
        self.edges      = defaultdict(set)
        self.in_degree  = defaultdict(int)
        self.out_degree = defaultdict(int)

    def _encode(self, s: str) -> int:
        base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        code = 0
        for b in s:
            code = (code << 2) | base_map[b]
        return code

    def _decode(self, code: int) -> str:
        inv = ['A', 'C', 'G', 'T']
        s = []
        for _ in range(self.k):
            s.append(inv[code & 3])
            code >>= 2
        return ''.join(reversed(s))

    def build(self, reads: List[str]):
        for read in reads:
            for i in range(len(read) - self.k):
                u = self._encode(read[i : i+self.k])
                v = self._encode(read[i+1 : i+self.k+1])
                self.edges[u].add(v)
                self.out_degree[u] += 1
                self.in_degree[v]  += 1

    def remove_dead_ends(self, max_depth: int = 5):
        all_nodes = set(self.edges) | set(self.in_degree) | set(self.out_degree)

        pred: defaultdict[int, Set[int]] = defaultdict(set)
        for u, succs in self.edges.items():
            for v in succs:
                pred[v].add(u)

        to_remove: Set[int] = set()
        depth: dict[int,int]  = {}

        queue = deque()
        for n in all_nodes:
            if self.out_degree.get(n, 0) == 0:
                queue.append(n)
                depth[n] = 0

        while queue:
            curr = queue.popleft()
            d    = depth[curr]
            if d >= max_depth:
                continue
            for p in pred[curr]:
                if self.out_degree.get(p, 0) == 1:
                    if p not in depth:
                        depth[p] = d + 1
                        queue.append(p)
                    to_remove.add(p)

        for n in to_remove:
            for s in list(self.edges.get(n, ())):
                self.in_degree[s] -= 1
            for p in pred.get(n, ()):
                self.out_degree[p] -= 1
                self.edges[p].discard(n)

            self.edges.pop(n,   None)
            self.in_degree.pop(n, None)
            self.out_degree.pop(n, None)

        del pred

    def remove_bubbles(self, max_depth: int = 100):
        pred = defaultdict(set)
        for u, succs in self.edges.items():
            for v in succs:
                pred[v].add(u)

        visited_pairs = set()
        to_remove = set()

        for src, succs in list(self.edges.items()):
            if len(succs) < 2:
                continue

            seen = {}
            paths = []
            queue = deque([(src, [src])])

            while queue:
                node, path = queue.popleft()
                if len(path) > max_depth:
                    continue
                for nb in self.edges.get(node, ()):
                    if nb in path:
                        continue
                    new_path = path + [nb]
                    if nb in seen:
                        paths.append((seen[nb], new_path))
                    else:
                        seen[nb] = new_path
                        queue.append((nb, new_path))

            for p1, p2 in paths:
                key = (tuple(p1), tuple(p2))
                if key in visited_pairs:
                    continue
                visited_pairs.add(key)
                visited_pairs.add((tuple(p2), tuple(p1)))
                weak = p1 if len(p1) < len(p2) else p2
                to_remove.update(weak)

        for n in to_remove:
            for s in list(self.edges.get(n, ())):
                self.in_degree[s] -= 1
            for p in pred.get(n, ()):
                if p in self.edges:
                    self.out_degree[p] -= 1
                    self.edges[p].discard(n)
            self.edges.pop(n, None)
            self.in_degree.pop(n, None)
            self.out_degree.pop(n, None)

        del pred
    
    def edge_count(self):
        return sum(len(succs) for succs in self.edges.values())


