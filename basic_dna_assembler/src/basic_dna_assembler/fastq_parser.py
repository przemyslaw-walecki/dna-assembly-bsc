from typing import List

class FastqParser:
    def __init__(self, file1: str, file2: str):
        self.file1 = file1
        self.file2 = file2

    def load_reads(self) -> List[str]:
        reads = []
        for file in [self.file1, self.file2]:
            with open(file) as f:
                while True:
                    header = f.readline()
                    if not header:
                        break
                    seq = f.readline().strip()
                    f.readline()
                    f.readline()
                    reads.append(seq)
        return reads
