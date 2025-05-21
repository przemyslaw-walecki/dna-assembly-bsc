from collections import defaultdict, deque
from typing import List


class KmerGraph:
    def __init__(self, k: int):
        """
        de Bruijn graph of k-mers with edge-coverage filtering.
        """
        self.k = k
        self.edges = defaultdict(set)  # u -> {v,…}
        self.in_degree = defaultdict(int)  # v -> count of incoming unique edges
        self.out_degree = defaultdict(int)  # u -> count of outgoing unique edges
        self.edge_counts = defaultdict(int)  # (u,v) -> no. of times seen

    def _encode(self, s: str) -> int:
        base = {"A": 0, "C": 1, "G": 2, "T": 3}
        c = 0
        for b in s:
            c = (c << 2) | base[b]
        return c

    def _decode(self, c: int) -> str:
        inv = ["A", "C", "G", "T"]
        s = []
        for _ in range(self.k):
            s.append(inv[c & 3])
            c >>= 2
        return "".join(reversed(s))

    def build(self, reads: List[str]):
        """
        Count every k→k+1 transition, record its coverage,
        then add it to the graph only once (unique edges).
        """
        for r in reads:
            for i in range(len(r) - self.k):
                u = self._encode(r[i : i + self.k])
                v = self._encode(r[i + 1 : i + self.k + 1])
                self.edge_counts[(u, v)] += 1
                if v not in self.edges[u]:
                    self.edges[u].add(v)
                    self.out_degree[u] += 1
                    self.in_degree[v] += 1

    def filter_low_coverage(self, threshold: int = 5):
        """
        Remove any edge seen fewer than `threshold` times.
        Adjusts in/out-degrees accordingly.
        """
        for (u, v), cnt in list(self.edge_counts.items()):
            if cnt < threshold and v in self.edges[u]:
                self.edges[u].remove(v)
                self.out_degree[u] -= 1
                self.in_degree[v] -= 1

    def remove_dead_ends(self, max_depth=None):
        if max_depth is None:
            max_depth = self.k

        for depth in range(1, max_depth + 1):

            pred = defaultdict(set)
            for u, vs in self.edges.items():
                for v in vs:
                    pred[v].add(u)

            to_remove = set()
            tips = [n for n, od in self.out_degree.items() if od == 0]
            for tip in tips:
                chain = [tip]
                curr = tip
                for _ in range(depth):
                    preds = [
                        p
                        for p in pred[curr]
                        if self.in_degree[p] == 1 and self.out_degree[p] == 1
                    ]
                    if len(preds) != 1:
                        break
                    curr = preds[0]
                    chain.append(curr)
                if len(chain) - 1 <= depth:
                    to_remove.update(chain)

            if not to_remove:
                break

            for n in to_remove:

                for v in list(self.edges.get(n, ())):
                    self.in_degree[v] -= 1

                for p in list(pred.get(n, ())):
                    if n in self.edges[p]:
                        self.out_degree[p] -= 1
                        self.edges[p].remove(n)
                self.edges.pop(n, None)
                self.in_degree.pop(n, None)
                self.out_degree.pop(n, None)

    def remove_bubbles(self, max_depth: int = 100):
        """
        Bubble removal: find pairs of distinct paths from the same source,
        drop the shorter one (up to max_depth).
        """
        pred = defaultdict(set)
        for u, vs in self.edges.items():
            for v in vs:
                pred[v].add(u)

        to_drop = set()
        seen_pairs = set()

        for src, vs in list(self.edges.items()):
            if len(vs) < 2:
                continue
            seen = {}
            q = deque([(src, [src])])
            while q:
                n, path = q.popleft()
                if len(path) > max_depth:
                    continue
                for nxt in self.edges.get(n, ()):
                    if nxt in path:
                        continue
                    newp = path + [nxt]
                    if nxt in seen:
                        p2 = seen[nxt]
                        key = (tuple(path), tuple(p2))
                        if key not in seen_pairs:
                            seen_pairs.add(key)
                            seen_pairs.add((tuple(p2), tuple(path)))

                            drop = path if len(path) < len(p2) else p2
                            to_drop.update(drop)
                    else:
                        seen[nxt] = newp
                        q.append((nxt, newp))

        for n in to_drop:
            for v in list(self.edges.get(n, ())):
                self.in_degree[v] -= 1
            for p in pred.get(n, ()):
                if n in self.edges[p]:
                    self.out_degree[p] -= 1
                    self.edges[p].discard(n)
            self.edges.pop(n, None)
            self.in_degree.pop(n, None)
            self.out_degree.pop(n, None)

    def edge_count(self) -> int:
        return sum(len(vs) for vs in self.edges.values())
