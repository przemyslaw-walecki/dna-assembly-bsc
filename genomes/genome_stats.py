import gzip
import json
from pathlib import Path
from typing import Dict, Any, Tuple, List

# ======================= KONFIGURACJA ==========================
# UZUPEŁNIJ TYLKO TO
CONFIG: Dict[str, Dict[str, Any]] = {
    "Ecoli": {
        "genome_fasta": "./Escherichia_coli_str_K-12_substr_MG1655.fasta",
        "coverages": {
            "10x": {
                "reads": [
                    "./ecoli_clean_cov10_100_300_1.fq",
                    "./ecoli_clean_cov10_100_300_2.fq",
                ],
                "gfa": "../rust_bubble_labelling/ecoli_k21_cov10_graph.gfa",
                "bubbles": "../bubble_resolution_gnn/lower_cov_data/ecoli_dataset_k21_cov10_ratio_lt_3.5.jsonl",
            },
            "20x": {
                "reads": [
                    "./ecoli_clean_cov20_100_300_1.fq",
                    "./ecoli_clean_cov20_100_300_2.fq",
                ],
                "gfa": "../rust_bubble_labelling/ecoli_k21_cov20_graph.gfa",
                "bubbles": "../bubble_resolution_gnn/lower_cov_data/ecoli_dataset_k21_cov20_ratio_lt_3.5.jsonl",
            },
            "30x": {
                "reads": [
                    "./ecoli_clean_cov30_100_300_1.fq",
                    "./ecoli_clean_cov30_100_300_2.fq",
                ],
                "gfa": "../rust_bubble_labelling/ecoli_k21_cov30_graph.gfa",
                "bubbles": "../bubble_resolution_gnn/lower_cov_data/ecoli_dataset_k21_cov30_ratio_lt_3.5.jsonl",
            },
        },
    },
    "Klebsiella": {
        "genome_fasta": "Klebsiella_pneumoniae_subsp_pneumoniae_NTUH-K2044_DNA.fasta",
        "coverages": {
            "10x": {
                "reads": [
                    "./klebsiella_clean_cov10_100_300_1.fq",
                    "./klebsiella_clean_cov10_100_300_2.fq",
                ],
                "gfa": "../rust_bubble_labelling/klebsiella_k21_cov10_graph.gfa",
                "bubbles": "../bubble_resolution_gnn/lower_cov_data/klebsiella_dataset_k21_cov10_ratio_lt_3.5.jsonl",
            },
            "20x": {
                "reads": [                    
                    "./klebsiella_clean_cov20_100_300_1.fq",
                    "./klebsiella_clean_cov20_100_300_2.fq",
                ],
                "gfa": "../rust_bubble_labelling/klebsiella_k21_cov20_graph.gfa",
                "bubbles": "../bubble_resolution_gnn/lower_cov_data/klebsiella_dataset_k21_cov20_ratio_lt_3.5.jsonl",
            },
            "30x": {
                "reads": [
                    "./klebsiella_clean_cov30_100_300_1.fq",
                    "./klebsiella_clean_cov30_100_300_2.fq",
                ],
                "gfa": "../rust_bubble_labelling/klebsiella_k21_cov30_graph.gfa",
                "bubbles": "../bubble_resolution_gnn/lower_cov_data/klebsiella_dataset_k21_cov30_ratio_lt_3.5.jsonl",
            },
        },
    },
    "Proteus": {
        "genome_fasta": "./Proteus-mirabilis.fasta",
        "coverages": {
            "10x": {
                "reads": [
                    "./proteus_clean_cov10_100_300_1.fq",
                    "./proteus_clean_cov10_100_300_2.fq",
                ],
                "gfa": "../rust_bubble_labelling/proteus_k21_cov10_graph.gfa",
                "bubbles": "../bubble_resolution_gnn/lower_cov_data/proteus_dataset_k21_cov10_ratio_lt_3.5.jsonl",
            },
            "20x": {
                "reads": [
                    "./proteus_clean_cov20_100_300_1.fq",
                    "./proteus_clean_cov20_100_300_2.fq",
                ],
                "gfa": "../rust_bubble_labelling/proteus_k21_cov20_graph.gfa",
                "bubbles": "../bubble_resolution_gnn/lower_cov_data/proteus_dataset_k21_cov20_ratio_lt_3.5.jsonl",
            },
            "30x": {
                "reads": [
                    "./proteus_clean_cov30_100_300_1.fq",
                    "./proteus_clean_cov30_100_300_2.fq",
                ],
                "gfa": "../rust_bubble_labelling/proteus_k21_cov30_graph.gfa",
                "bubbles": "../bubble_resolution_gnn/lower_cov_data/proteus_dataset_k21_cov30_ratio_lt_3.5.jsonl",
            },
        },
    },
}
# ===============================================================


def open_maybe_gzip(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def genome_size_fasta(path: Path) -> int:
    size = 0
    with open_maybe_gzip(path, "rt") as f:
        for line in f:
            if not line:
                continue
            if line[0] == ">":
                continue
            size += len(line.strip())
    return size


def count_reads_fastq(path: Path) -> int:
    # zakładamy standardowy FASTQ: 4 linie na odczyt
    lines = 0
    with open_maybe_gzip(path, "rt") as f:
        for _ in f:
            lines += 1
    return lines // 4


def count_nodes_edges_gfa(path: Path) -> Tuple[int, int]:
    nodes = 0
    edges = 0
    with open_maybe_gzip(path, "rt") as f:
        for line in f:
            if not line:
                continue
            c = line[0]
            if c == "S":  # segment
                nodes += 1
            elif c == "L":  # link
                edges += 1
    return nodes, edges


def count_bubbles(path: Path) -> int:
    # JSONL: 1 obiekt na linię (Twoje wcześniejsze ecoli_dataset_...jsonl)
    if path.suffix == ".jsonl":
        cnt = 0
        with open_maybe_gzip(path, "rt") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                json.loads(line)  # walidacja
                cnt += 1
        return cnt

    # JSON: słownik/lub lista
    with open_maybe_gzip(path, "rt") as f:
        obj = json.load(f)

    # przypadek: lista superbabli
    if isinstance(obj, list):
        return len(obj)

    # przypadek: słownik chain_id -> {"bubbles":[...]}
    if isinstance(obj, dict):
        total = 0
        for _, chain in obj.items():
            bubbles = chain.get("bubbles", [])
            # licz tylko superbable (u Ciebie wszystkie są "super", ale filtr jest bezpieczniejszy)
            total += sum(1 for b in bubbles if b.get("type") == "super")
        return total

    raise ValueError(f"Nieobsługiwany format JSON w pliku: {path}")
def main():
    # policz rozmiary genomów raz
    genome_sizes: Dict[str, int] = {}
    for genome_name, cfg in CONFIG.items():
        fasta_path = Path(cfg["genome_fasta"])
        genome_sizes[genome_name] = genome_size_fasta(fasta_path)

    # przygotuj tabelę
    rows: List[List[str]] = []
    header = [
        "genome",
        "genome_size_bp",
        "coverage",
        "num_reads",
        "graph_nodes",
        "graph_edges",
        "num_superbubbles",
    ]
    rows.append(header)

    for genome_name, cfg in CONFIG.items():
        for cov, cov_cfg in cfg["coverages"].items():
            # odczyty
            total_reads = 0
            for rpath in cov_cfg["reads"]:
                total_reads += count_reads_fastq(Path(rpath))

            # graf
            nodes, edges = count_nodes_edges_gfa(Path(cov_cfg["gfa"]))

            # superbable
            superbubbles = count_bubbles(Path(cov_cfg["bubbles"]))

            rows.append(
                [
                    genome_name,
                    str(genome_sizes[genome_name]),
                    cov,
                    str(total_reads),
                    str(nodes),
                    str(edges),
                    str(superbubbles),
                ]
            )

    for row in rows:
        print("\t".join(row))


if __name__ == "__main__":
    main()
