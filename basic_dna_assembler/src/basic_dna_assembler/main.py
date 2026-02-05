"""
Główny moduł uruchamiający podstawowy assembler DNA w Pythonie.

Realizuje pipeline:
1. Wczytanie odczytów z dwóch plików FASTQ.
2. Zbudowanie grafu de Bruijna dla zadanej długości k-merów.
3. Filtrowanie krawędzi o niskim pokryciu.
4. Usuwanie krótkich zakończeń (tips).
5. Usuwanie prostych baniek.
6. Składanie kontigów.
7. Zapis wyników do pliku FASTA w katalogu `./data`.
"""

import argparse
import os
from src.basic_dna_assembler.fastq_parser import FastqParser
from src.basic_dna_assembler.kmer_graph import KmerGraph
from src.basic_dna_assembler.assembler import Assembler


def main():
    """
    Uruchamia assembler na podstawie argumentów wiersza poleceń.

    Obsługiwane opcje:
        -1 / --reads1  : plik FASTQ z pierwszymi końcami odczytów,
        -2 / --reads2  : plik FASTQ z drugimi końcami,
        -k / --kmer    : długość k-mera,
        -t / --threshold : minimalne pokrycie krawędzi,
        -d / --depth   : maksymalna długość usuwanych zakończeń,
        -o / --output  : nazwa pliku wynikowego FASTA.

    Zwraca:
        Plik FASTA z kontigami zapisany w `./data/`.
    """
    p = argparse.ArgumentParser(
        description="Assembler genomowy oparty na grafie de Bruijna."
    )
    p.add_argument("-1", "--reads1", required=True)
    p.add_argument("-2", "--reads2", required=True)
    p.add_argument("-k", "--kmer", type=int, default=55)
    p.add_argument("-t", "--threshold", type=int, default=5)
    p.add_argument("-o", "--output", default="assembled.fa")
    p.add_argument("-d", "--depth", type=int, default=5)

    args = p.parse_args()

    reads1 = FastqParser(args.reads1).load_reads()
    reads2 = FastqParser(args.reads2).load_reads()
    reads = reads1 + reads2

    graph = KmerGraph(args.kmer)
    graph.build(reads)
    graph.filter_low_coverage(threshold=args.threshold)
    graph.remove_dead_ends(max_depth=args.depth)
    graph.remove_bubbles(max_depth=20)

    contigs = Assembler(graph).assemble_contigs()

    out_path = os.path.join("data", args.output)
    with open(out_path, "w") as f:
        for i, c in enumerate(contigs, start=1):
            f.write(f">contig_{i}\n{c}\n")

    print(f"Zapisano {len(contigs)} kontigów do {out_path}")


if __name__ == "__main__":
    main()
