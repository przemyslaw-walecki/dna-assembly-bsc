import argparse
import os
from src.basic_dna_assembler.fastq_parser import FastqParser
from src.basic_dna_assembler.kmer_graph import KmerGraph
from src.basic_dna_assembler.assembler import Assembler


def main():
    p = argparse.ArgumentParser(
        description="De Bruijn-based assembler using two FASTQ files"
    )
    p.add_argument("-1", "--reads1", required=True, help="First-end FASTQ file")
    p.add_argument("-2", "--reads2", required=True, help="Second-end FASTQ file")
    p.add_argument("-k", "--kmer", type=int, default=55, help="k-mer size")
    p.add_argument(
        "-t", "--threshold", type=int, default=5, help="Min edge coverage to keep"
    )
    p.add_argument(
        "-o",
        "--output",
        default="assembled.fa",
        help="Output FASTA filename (written into ./data/)",
    )
    p.add_argument(
        "-d", "--depth", type=int, default=5, help="Depth of dead ends to remove"
    )
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

    print(f"Wrote {len(contigs)} contigs to {out_path}")


if __name__ == "__main__":
    main()
