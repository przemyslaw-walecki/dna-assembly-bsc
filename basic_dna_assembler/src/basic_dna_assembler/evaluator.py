import subprocess
import os
import re

class Evaluator:
    def __init__(self, reference_path: str, contigs_path: str, output_sam="alignment.sam"):
        self.reference_path = reference_path
        self.contigs_path = contigs_path
        self.output_sam = output_sam

    def evaluate(self):
        print("Evaluating assembly")
        cmd = [
            "minimap2", "-a", self.reference_path, self.contigs_path
        ]
        try:
            with open(self.output_sam, "w") as sam_out:
                subprocess.run(cmd, stdout=sam_out, stderr=subprocess.PIPE, text=True, check=True)
        except FileNotFoundError:
            print("minimap2 not installed. Skipping evaluation.")
            return
        except subprocess.CalledProcessError as e:
            print(f"Error during minimap2 execution:\n{e.stderr}")
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

        print("Assembly evaluation summary:")
        print(f"  Total contigs aligned: {aligned}/{total}")
        if lengths:
            avg_len = sum(lengths) / len(lengths)
            print(f"  Avg aligned segment length: {avg_len:.1f} bp")
        if mismatches:
            avg_mismatch = sum(mismatches) / len(mismatches)
            print(f"  Avg mismatches per alignment: {avg_mismatch:.2f}")
        print(f"  Full SAM saved to: {self.output_sam}")
