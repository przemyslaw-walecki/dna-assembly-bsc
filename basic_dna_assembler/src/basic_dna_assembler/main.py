from src.basic_dna_assembler.assembler import Assembler
from src.basic_dna_assembler.evaluator import Evaluator
from src.basic_dna_assembler.fastq_parser import FastqParser
from src.basic_dna_assembler.kmer_graph import KmerGraph

def main():
    reads_1 = "./data/ecoli_reads_100_400_1.fq"
    reads_2 = "./data/ecoli_reads_100_400_2.fq"
    reference = "./data/GCF_000007365.1_ASM736v1_genomic.fna"

    parser = FastqParser(reads_1, reads_2)
    reads = parser.load_reads()


    graph = KmerGraph(k=51)
    graph.build(reads)
    print(F"Built graph. Edges: {graph.edge_count()}.")
    graph.remove_dead_ends(20)
    print(f"Removed dead branches. Edges after pruning: {graph.edge_count()}.")
    branching = sum(1 for succs in graph.edges.values() if len(succs) > 1)
    print(f" {branching} nodes with out_degree > 1 (branching)")
    #graph.remove_bubbles(10)
    #print(f"Removed bubbles. Edges for asssembly: {graph.edge_count()}.")

    assembler = Assembler(graph)
    contigs = assembler.assemble_contigs()

    with open("assembled_contigs.fasta", "w") as f:
        for i, contig in enumerate(contigs):
            f.write(f">contig_{i}\n{contig}\n")

    evaluator = Evaluator(reference, "assembled_contigs.fasta")
    evaluator.evaluate()

if __name__ == "__main__":
    main()
