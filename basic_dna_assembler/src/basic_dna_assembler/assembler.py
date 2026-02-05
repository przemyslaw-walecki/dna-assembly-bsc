"""
Moduł odpowiedzialny za składanie kontigów na podstawie grafu de Bruijna.

Zawiera klasę :class:`Assembler`, która przetwarza oczyszczony graf k-merów
i odtwarza sekwencje DNA poprzez odnajdywanie maksymalnych ścieżek
niebranchujących (segmentów 1→1). Każda taka ścieżka jest dekodowana do
ciągu nukleotydów i zapisywana jako kontig.
"""


class Assembler:
    """
    Prosty assembler kontigów działający na grafie de Bruijna.

    Algorytm wyszukuje wszystkie maksymalne ścieżki, w których każdy węzeł
    wewnętrzny ma dokładnie jeden poprzednik i jednego następnika. Pierwszy
    k-mer w ścieżce jest dekodowany w całości, a każdy kolejny k-mer
    wydłuża sekwencję o jedną zasadę.
    """

    def __init__(self, graph):
        """
        Inicjalizuje assembler.

        Argumenty:
            graph: Oczyszczony graf de Bruijna zawierający:
                - `edges`: mapowanie węzeł → zbiór następców,
                - `in_degree`: liczba krawędzi wchodzących,
                - `out_degree`: liczba krawędzi wychodzących,
                - `_decode(int) -> str`: metodę dekodującą k-mery.
        """
        self.graph = graph

    def assemble_contigs(self):
        """
        Składa kontigi poprzez przechodzenie maksymalnych ścieżek niebranchujących.

        Każda ścieżka zaczyna się w węźle, który nie ma stopnia wejścia = 1
        lub stopnia wyjścia = 1. Następnie ścieżka jest wydłużana, dopóki
        kolejne węzły spełniają warunek 1→1.

        Zwraca:
            list[str]: Lista zrekonstruowanych kontigów.
        """
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
