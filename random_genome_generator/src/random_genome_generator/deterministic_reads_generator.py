"""
generate_reads.py: Deterministic sliding-window read generator

Reads a FASTA file and generates all possible reads of a specified length
in FASTQ format with a fixed quality score ('I').

Usage:
    python generate_reads.py <source_fasta> <read_length> [--output OUTPUT]

Example:
    python generate_reads.py random_50.fna 20 --output reads50.fq
"""
import argparse
from pathlib import Path

def parse_fasta(path):
    """Simple FASTA parser yielding (header, sequence) tuples."""
    header = None
    seq_lines = []
    with open(path, 'r') as f:
        for line in f:
            line = line.rstrip()  # remove newline
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    yield header, ''.join(seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            yield header, ''.join(seq_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Generate deterministic sliding-window reads from a FASTA.")
    parser.add_argument(
        'source',
        help='Input FASTA file containing one or more reference sequences')
    parser.add_argument(
        'read_length',
        type=int,
        help='Length of each read to generate (must be <= reference length)')
    parser.add_argument(
        '--output',
        '-o',
        help='Output FASTQ file (default: <source_stem>_reads_<read_length>.fq)')
    args = parser.parse_args()

    read_len = args.read_length
    src = Path(args.source)
    out_path = Path(args.output) if args.output else src.with_name(f"{src.stem}_reads_{read_len}.fq")

    total_reads = 0
    with open(out_path, 'w') as out_f:
        for header, seq in parse_fasta(src):
            n = len(seq)
            if read_len > n:
                raise ValueError(f"Read length {read_len} is longer than sequence {header} length {n}")
            for i in range(n - read_len + 1):
                read_seq = seq[i:i + read_len]
                quals = 'I' * read_len
                out_f.write(f"@{header}_{i}\n")
                out_f.write(f"{read_seq}\n")