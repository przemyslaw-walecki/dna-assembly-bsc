"""
Moduł narzędziowy do porównywania wyników składania genomów.

Umożliwia:
- wyrównanie kontigów do genomu referencyjnego za pomocą programu minimap2,
- obliczenie statystyk dopasowania,
- zapisanie wyników do pliku CSV,
- generowanie podsumowania w formacie Markdown.

Zawiera także interfejs CLI, który może analizować wiele assemblerów naraz
przy założeniu, że kontigi zostały zapisane jako `assembled_{nazwa}.fa`.
"""

import argparse
import os
import csv
import textwrap
import subprocess
import re


def evaluate_assembly(name, reference, contigs_path, writer, ref_label):
    """
    Wyrównuje kontigi do genomu referencyjnego i zapisuje statystyki.

    Funkcja uruchamia minimap2 w trybie SAM, analizuje wynik i oblicza:
    - liczbę wszystkich rekordów,
    - liczbę rekordów dopasowanych,
    - procent dopasowanych,
    - średnią długość segmentów typu „M” w CIGAR,
    - średnią liczbę błędów (znacznik NM).

    Wyniki są zapisywane jako jeden wiersz w pliku CSV.

    Argumenty:
        name (str): Nazwa assemblera.
        reference (str): Ścieżka do genomu referencyjnego.
        contigs_path (str): Ścieżka do pliku FASTA z kontigami.
        writer (csv.DictWriter): Obiekt zapisujący wiersz danych.
        ref_label (str): Etykieta opisująca genom referencyjny.
    """
    sam_file = os.path.abspath(f"{name}.sam")
    try:
        proc = subprocess.run(
            ["minimap2", "-a", reference, contigs_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    except FileNotFoundError:
        print("Błąd: minimap2 nie jest zainstalowany.")
        return
    except subprocess.CalledProcessError as e:
        print(f"Błąd podczas uruchamiania minimap2 dla {name}:\n{e.stderr}")
        return

    with open(sam_file, "w") as f:
        f.write(proc.stdout)

    total = aligned = 0
    lengths = []
    mismatches = []

    for line in proc.stdout.splitlines():
        if line.startswith("@"):
            continue
        total += 1
        fields = line.split("\t")
        flag = int(fields[1])
        if flag & 4:
            continue
        aligned += 1

        for m in re.findall(r"(\d+)M", fields[5]):
            lengths.append(int(m))
        for tag in fields[11:]:
            if tag.startswith("NM:i:"):
                mismatches.append(int(tag.split(":")[-1]))

    avg_len = sum(lengths) / len(lengths) if lengths else 0.0
    avg_mm = sum(mismatches) / len(mismatches) if mismatches else 0.0
    pct = aligned / total * 100 if total else 0.0

    writer.writerow(
        {
            "reference": ref_label,
            "assembler": name,
            "total": total,
            "aligned": aligned,
            "pct_aligned": f"{pct:.1f}",
            "avg_len": f"{avg_len:.1f}",
            "avg_mm": f"{avg_mm:.2f}",
        }
    )

    print(
        f"-> {name}: {aligned}/{total} dopasowanych "
        f"(śr. długość={avg_len:.1f}, śr. błędy={avg_mm:.2f})"
    )


def append_markdown_table(csv_path, md_log, ref_label):
    """
    Dodaje do pliku Markdown tabelę podsumowującą wyniki z pliku CSV.

    Argumenty:
        csv_path (str): Ścieżka do pliku CSV z wynikami.
        md_log (str): Ścieżka do pliku Markdown.
        ref_label (str): Opis analizowanego genomu referencyjnego.
    """
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))

    headers = [
        "Referencja",
        "Assembler",
        "Łącznie",
        "Dopasowane",
        "% Dopas.",
        "Śr. długość",
        "Śr. błędy",
    ]

    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]

    for r in rows:
        lines.append(
            f"| {r['reference']} | {r['assembler']} | {r['total']} "
            f"| {r['aligned']} | {r['pct_aligned']} | {r['avg_len']} | {r['avg_mm']} |"
        )

    with open(md_log, "a") as log:
        log.write(
            textwrap.dedent(
                f"""

        ### Wyniki porównania montażu

        **Genom referencyjny:** `{ref_label}`

        """
            )
        )
        log.write("\n".join(lines) + "\n")


def main():
    """
    Główny punkt wejścia CLI do porównywania wyników montażu genomów.

    Skrypt:
    1. Sprawdza istnienie genomu referencyjnego.
    2. Dla każdego assemblera wyszukuje plik `assembled_{nazwa}.fa`.
    3. Wykonuje wyrównanie minimap2.
    4. Zapisuje wyniki do CSV.
    5. Dopisuje podsumowanie w formacie Markdown.

    Obsługiwane opcje:
        --data-dir  : katalog z kontigami,
        --reference : genom referencyjny,
        --assemblers: lista nazw assemblerów,
        --output    : wynikowy plik CSV,
        --md-log    : log w formacie Markdown.
    """
    parser = argparse.ArgumentParser(
        description="Ewaluacja i porównanie gotowych plików kontigów."
    )
    parser.add_argument("--data-dir", "-d", default="./data")
    parser.add_argument("--reference", "-r", required=True)
    parser.add_argument("--assemblers", "-a", nargs="+", default=["dnaasm"])
    parser.add_argument("--output", "-o", default="assembly_comparison.csv")
    parser.add_argument("--md-log", "-m", default="comparison.log.md")

    args = parser.parse_args()

    reference = os.path.abspath(args.reference)
    if not os.path.exists(reference):
        parser.error(f"Nie znaleziono referencji: {reference}")
    ref_label = "/" + os.path.basename(reference)

    with open(args.output, "w", newline="") as csvfile:
        fieldnames = [
            "reference",
            "assembler",
            "total",
            "aligned",
            "pct_aligned",
            "avg_len",
            "avg_mm",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for name in args.assemblers:
            contigs = os.path.join(args.data_dir, f"assembled_{name}.fa")
            if not os.path.exists(contigs):
                print(f"→ Pomijam {name}: brak pliku {contigs}")
                continue
            evaluate_assembly(name, reference, contigs, writer, ref_label)

    append_markdown_table(args.output, args.md_log, ref_label)
    print(f"\nWyniki zapisane do {args.output}")
    print(f"Podsumowanie dopisano do {args.md_log}")


if __name__ == "__main__":
    main()
