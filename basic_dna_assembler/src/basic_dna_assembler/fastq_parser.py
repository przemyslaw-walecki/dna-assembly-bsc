"""
Moduł do prostego wczytywania sekwencji z plików FASTQ.

Parser uwzględnia wyłącznie sekwencje z pierwszego pliku pary FASTQ,
co pozwala uniknąć duplikacji wynikającej z odwróconych komplementów
w sztucznie generowanych genomach i jest wystarczające na potrzeby
podstawowego assemblera.
"""

from typing import List


class FastqParser:
    """
    Minimalistyczny parser FASTQ do odczytu sekwencji nukleotydowych.

    Odczytuje standardowe czteroliniowe rekordy FASTQ i zwraca listę
    sekwencji DNA bez nagłówków i bez linii jakości.
    """

    def __init__(self, file1: str):
        """
        Inicjalizuje parser.

        Argumenty:
            file1 (str): Ścieżka do pliku FASTQ.
        """
        self.file1 = file1

    def load_reads(self) -> List[str]:
        """
        Wczytuje wszystkie sekwencje z pliku FASTQ.

        Zwraca:
            list[str]: Lista odczytanych sekwencji DNA.
        """
        reads = []
        with open(self.file1) as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()  # linia "+"
                f.readline()  # linia jakości
                reads.append(seq)
        return reads
