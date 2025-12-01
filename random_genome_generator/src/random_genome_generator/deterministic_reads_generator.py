"""
Generator deterministycznych odczytów sliding-window z pliku FASTA.

Moduł wczytuje sekwencje referencyjne w formacie FASTA i generuje wszystkie
możliwe odczyty o zadanej długości, przesuwając okno o jeden nukleotyd.
Wynik zapisywany jest w formacie FASTQ z jednolitym ciągiem jakości ('I').

Przykładowe użycie:
    python deterministic_reads_generator.py source.fna 20 --output reads.fq

Zastosowanie:
    - testowanie assemblerów,
    - walidacja pipeline’ów montażu,
    - generowanie syntetycznych danych pozbawionych losowości.
"""

import argparse
from pathlib import Path


def parse_fasta(path):
    """
    Prosty parser FASTA zwracający pary (nagłówek, sekwencja).

    Działa strumieniowo — nie wczytuje całego pliku do pamięci.

    Argumenty:
        path (str | Path): Ścieżka do pliku FASTA.

    Zwraca:
        generator[tuple[str, str]]: Para (header, sequence).
    """
    header = None
    seq_lines = []
    with open(path, 'r') as f:
        for line in f:
            line = line.rstrip()
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
    """
    Główna funkcja CLI generująca deterministyczne odczyty FASTQ.

    Działanie:
        1. Wczytuje sekwencje z pliku FASTA.
        2. Dla każdej sekwencji generuje wszystkie możliwe odczyty
           o zadanej długości (okno przesuwane o 1).
        3. Zapisuje odczyty w formacie FASTQ z jakością 'I'.

    Argumenty wejściowe:
        source (plik FASTA),
        read_length (długość odczytu),
        --output / -o (ścieżka wyjściowa FASTQ).

    Wynik:
        Plik FASTQ zawierający deterministyczne odczyty sliding-window.
    """
    parser = argparse.ArgumentParser(
        description="Generuj deterministyczne odczyty sliding-window z FASTA."
    )
    parser.add_argument(
        'source',
        help='Plik FASTA zawierający jedną lub więcej sekwencji referencyjnych'
    )
    parser.add_argument(
        'read_length',
        type=int,
        help='Długość pojedynczego odczytu (musi być ≤ długości sekwencji)'
    )
    parser.add_argument(
        '--output',
        '-o',
        help='Plik wyjściowy FASTQ (domyślnie <nazwa>_reads_<len>.fq)'
    )
    args = parser.parse_args()

    read_len = args.read_length
    src = Path(args.source)
    out_path = Path(args.output) if args.output else src.with_name(
        f"{src.stem}_reads_{read_len}.fq"
    )

    with open(out_path, 'w') as out_f:
        for header, seq in parse_fasta(src):
            n = len(seq)
            if read_len > n:
                raise ValueError(
                    f"Długość odczytu {read_len} przekracza długość sekwencji {header} ({n})"
                )
            for i in range(n - read_len + 1):
                read_seq = seq[i:i + read_len]
                quals = 'I' * read_len
                out_f.write(f"@{header}_{i}\n")
                out_f.write(f"{read_seq}\n")
                out_f.write("+\n")
                out_f.write(f"{quals}\n")
