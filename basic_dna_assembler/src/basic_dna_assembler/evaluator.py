"""
Moduł do podstawowej oceny jakości montażu genomu.

Pozwala na:
- uruchomienie programu minimap2 dla pary: plik kontigów + genom referencyjny,
- zapisanie pełnego pliku SAM,
- obliczenie podstawowych statystyk wyrównania,
- wyświetlenie krótkiego podsumowania na standardowym wyjściu.

Moduł stanowi uproszczony wariant narzędzia porównawczego, przydatny
do szybkiej walidacji pojedynczego montażu.
"""

import subprocess
import re


class Evaluator:
    """
    Klasa służąca do oceny jakości montażu poprzez analizę wyrównania SAM.

    Uruchamia minimap2, generuje plik SAM, a następnie oblicza m.in. liczbę
    dopasowanych rekordów, średnią długość segmentu wyrównanego i średnią
    liczbę błędów NM.
    """

    def __init__(self, reference_path: str, contigs_path: str, output_sam="alignment.sam"):
        """
        Inicjalizuje obiekt evaluator.

        Argumenty:
            reference_path (str): Ścieżka do genomu referencyjnego.
            contigs_path (str): Ścieżka do pliku FASTA z kontigami.
            output_sam (str): Nazwa pliku wynikowego SAM.
        """
        self.reference_path = reference_path
        self.contigs_path = contigs_path
        self.output_sam = output_sam

    def evaluate(self):
        """
        Wykonuje wyrównanie minimap2 i oblicza statystyki.

        Działanie:
            1. Uruchamia `minimap2 -a`.
            2. Zapisuje pełen wynik wyrównania w pliku SAM.
            3. Analizuje rekordy wyrównane i oblicza podstawowe statystyki.
            4. Wyświetla podsumowanie.

        Jeśli minimap2 nie jest dostępny, metoda kończy działanie z komunikatem.
        """
        print("Rozpoczynam ocenę montażu.")
        cmd = ["minimap2", "-a", self.reference_path, self.contigs_path]
        try:
            with open(self.output_sam, "w") as sam_out:
                subprocess.run(
                    cmd, stdout=sam_out, stderr=subprocess.PIPE, text=True, check=True
                )
        except FileNotFoundError:
            print("minimap2 nie jest zainstalowany — pomijam ocenę.")
            return
        except subprocess.CalledProcessError as e:
            print(f"Błąd podczas wykonywania minimap2:\n{e.stderr}")
            return

        aligned = 0
        total = 0
        lengths = []
        mismatches = []

        with open(self.output_sam, "r") as f:
            for line in f:
                if line.startswith("@"):
                    continue
                total += 1
                fields = line.strip().split("\t")
                flag = int(fields[1])
                if flag & 4:
                    continue
                aligned += 1

                cigar = fields[5]
                m = re.match(r"(\d+)M", cigar)
                if m:
                    lengths.append(int(m.group(1)))

                for field in fields[11:]:
                    if field.startswith("NM:i:"):
                        mismatches.append(int(field.split(":")[-1]))

        print("Podsumowanie oceny montażu:")
        print(f"  Dopasowane rekordy: {aligned}/{total}")
        if lengths:
            print(f"  Średnia długość dopasowania: {sum(lengths)/len(lengths):.1f} bp")
        if mismatches:
            print(f"  Średnia liczba błędów NM: {sum(mismatches)/len(mismatches):.2f}")
        print(f"  Zapisano pełny plik SAM: {self.output_sam}")
