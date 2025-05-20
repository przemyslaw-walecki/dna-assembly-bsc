import argparse
import random


def generate_random_genome(length):
    return ''.join(random.choices('ACGT', k=length))


def write_fasta(seq, header, output_path, wrap=60):
    with open(output_path, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), wrap):
            f.write(seq[i:i+wrap] + "\n")


def main():
    parser = argparse.ArgumentParser("generate-genome")
    parser.add_argument('--length', type=int, required=True,
                        help="Length of random genome to generate (bp)")
    parser.add_argument('--output', '-o', type=str, required=True,
                        help="Output FASTA file path")
    parser.add_argument('--wrap', type=int, default=60,
                        help="Line wrap length for FASTA (default: 60)")
    args = parser.parse_args()

    seq = generate_random_genome(args.length)
    header = f"random_genome_{args.length}bp"
    write_fasta(seq, header, args.output, wrap=args.wrap)
    print(f"Generated random genome ({args.length} bp) -> {args.output}")


if __name__ == '__main__':
    main()
