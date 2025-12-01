"""
Generator losowych genomów na potrzeby testów assemblerów i pipeline’ów NGS.

Moduł umożliwia:
- generowanie w pełni losowej sekwencji genomowej o zadanej długości,
- zapis wyniku w formacie FASTA z opcjonalnym łamaniem linii.

Zastosowania:
    - testowanie algorytmów montażu,
    - walidacja jakości odczytów,
    - generowanie syntetycznych benchmarków.
"""

import argparse
import random


def generate_random_genome(length):
    """
    Generuje losową sekwencję genomu o podanej długości.

    Argumenty:
        length (int): Liczba nukleotydów.

    Zwraca:
        str: Sekwencja złożona z liter A/C/G/T.
    """
    return ''.join(random.choices('ACGT', k=length))


def write_fasta(seq, header, output_path, wrap=60):
    """
    Zapisuje sekwencję do pliku FASTA.

    Argumenty:
        seq (str): Sekwencja genomowa.
        header (str): Nazwa sekwencji (bez znaku '>').
        output_path (str): Ścieżka do pliku wynikowego.
        wrap (int): Liczba znaków w linii (domyślnie 60).
    """
    with open(output_path, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), wrap):
            f.write(seq[i:i+wrap] + "\n")


def main():
    """
    Punkt wejścia CLI do generatora losowych genomów.

    Parametry:
        --length (int): Długość generowanej sekwencji.
        --output (-o): Plik wynikowy FASTA.
        --wrap (int): Długość linii w FASTA (domyślnie 60).

    Wynik:
        Zapisany plik FASTA z losową sekwencją genomową.
    """
    parser = argparse.ArgumentParser("generate-genome")
    parser.add_argument('--length', type=int, required=True,
                        help="Długość losowego genomu (bp)")
    parser.add_argument('--output', '-o', type=str, required=True,
                        help="Ścieżka pliku wyjściowego FASTA")
    parser.add_argument('--wrap', type=int, default=60,
                        help="Długość linii FASTA (domyślnie: 60)")
    args = parser.parse_args()

    seq = generate_random_genome(args.length)
    header = f"random_genome_{args.length}bp"
    write_fasta(seq, header, args.output, wrap=args.wrap)
    print(f"Wygenerowano genom losowy ({args.length} bp) → {args.output}")


if __name__ == '__main__':
    main()
