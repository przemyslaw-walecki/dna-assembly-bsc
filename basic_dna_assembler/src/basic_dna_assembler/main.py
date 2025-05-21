from src.basic_dna_assembler.assembler import Assembler
from src.basic_dna_assembler.fastq_parser import FastqParser
from src.basic_dna_assembler.kmer_graph import KmerGraph

def main():
    reads_1 = "../random_genome_generator/reads_100_400_1.fq"
    reads_2 = "../random_genome_generator/reads_100_400_2.fq"
    out_path = "./data/own_assembled_random.fa"

    parser = FastqParser(reads_1, reads_2)
    reads = parser.load_reads()

    graph = KmerGraph(k=51)
    graph.build(reads)
    print(F"Built graph. Edges: {graph.edge_count()}.")
    assembler = Assembler(graph)

    contigs = assembler.assemble_contigs()
    with open(out_path, "w") as f:
        for i, contig in enumerate(contigs):
            f.write(f">contig_{i}\n{contig}\n")

    print(f"Assembled {len(contigs)} contigs.")

if __name__ == "__main__":
    main()
