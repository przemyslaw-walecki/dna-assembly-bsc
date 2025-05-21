from typing import List


class FastqParser:
    def __init__(self, file1: str):
        """
        Only parse the first-end FASTQ to avoid
        reverse-complement duplication in random genomes.
        """
        self.file1 = file1

    def load_reads(self) -> List[str]:
        reads = []
        with open(self.file1) as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                reads.append(seq)
        return reads
