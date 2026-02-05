"""
Implementacja grafu de Bruijna opartego na k-merach.

Węzły grafu reprezentowane są jako liczby całkowite, kodowane w systemie
2-bitowym na nukleotyd. Moduł umożliwia:
- budowę grafu na podstawie zestawu odczytów,
- filtrowanie krawędzi o niskim pokryciu,
- usuwanie krótkich zakończeń (dead-endów),
- uproszczone usuwanie baniek (bubble popping),
- zliczanie krawędzi.

Graf jest stosowany w podstawowym assemblerze kontigów.
"""

from collections import defaultdict, deque
from typing import List


class KmerGraph:
    """
    Graf de Bruijna z k-merów reprezentowanych jako liczby całkowite.

    Zawiera podstawowe funkcje czyszczenia grafu oraz metodę budowy struktury
    na podstawie listy odczytów.
    """

    def __init__(self, k: int):
        """
        Inicjalizuje pusty graf de Bruijna.

        Argumenty:
            k (int): Długość k-merów.
        """
        self.k = k
        self.edges = defaultdict(set)
        self.in_degree = defaultdict(int)
        self.out_degree = defaultdict(int)
        self.edge_counts = defaultdict(int)

    def _encode(self, s: str) -> int:
        """
        Koduje k-mer do liczby całkowitej (2 bity na nukleotyd).

        Argumenty:
            s (str): K-mer o długości k.

        Zwraca:
            int: Liczbową reprezentację k-mera.
        """
        base = {"A": 0, "C": 1, "G": 2, "T": 3}
        c = 0
        for b in s:
            c = (c << 2) | base[b]
        return c

    def _decode(self, c: int) -> str:
        """
        Dekoduje liczbową reprezentację k-mera do formy tekstowej.

        Argumenty:
            c (int): Liczba zakodowana metodą `_encode`.

        Zwraca:
            str: K-mer w postaci tekstowej (A/C/G/T).
        """
        inv = ["A", "C", "G", "T"]
        s = []
        for _ in range(self.k):
            s.append(inv[c & 3])
            c >>= 2
        return "".join(reversed(s))

    def build(self, reads: List[str]):
        """
        Buduje graf de Bruijna na podstawie listy odczytów.

        Dla każdego odczytu wykrywane są wszystkie k→k+1 przejścia
        i zapisywane jako krawędzie oraz zliczane jako pokrycie.

        Argumenty:
            reads (list[str]): Lista sekwencji z FASTQ.
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
        Usuwa krawędzie o zbyt niskim pokryciu.

        Argumenty:
            threshold (int): Minimalne pokrycie wymagane do pozostawienia
                krawędzi. Domyślnie 5.
        """
        for (u, v), cnt in list(self.edge_counts.items()):
            if cnt < threshold and v in self.edges[u]:
                self.edges[u].remove(v)
                self.out_degree[u] -= 1
                self.in_degree[v] -= 1

    def remove_dead_ends(self, max_depth=None):
        """
        Usuwa krótkie zakończenia (tzw. tips) z grafu.

        Algorytm:
            - wyszukuje węzły o stopniu wyjścia 0,
            - cofa się po łańcuchach 1→1,
            - jeśli długość łańcucha nie przekracza max_depth, usuwa go.

        Argumenty:
            max_depth (int | None): Maksymalna długość usuwalnego zakończenia.
                Jeśli `None`, przyjmowana jest wartość `k`.
        """
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
        Usuwa proste bańki (bubble popping) z grafu.

        Bańką nazywamy rozgałęzienie, które ponownie łączy się w późniejszym
        węźle. Algorytm wykrywa dwie różne ścieżki dochodzące do tego samego
        węzła i usuwa tę krótszą.

        Argumenty:
            max_depth (int): Maksymalna analizowana długość ścieżek.
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
        """
        Zwraca łączną liczbę krawędzi w grafie.

        Zwraca:
            int: Liczbę krawędzi skierowanych.
        """
        return sum(len(vs) for vs in self.edges.values())
