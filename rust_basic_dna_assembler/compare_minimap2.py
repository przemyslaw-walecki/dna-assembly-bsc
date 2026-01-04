#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess
import re
import textwrap


def run_minimap2(reference, assembly_path):
    """
    Run minimap2 -a (SAM output) and return SAM text.
    """
    try:
        proc = subprocess.run(
            ["minimap2", "-a", reference, assembly_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        return proc.stdout
    except FileNotFoundError:
        raise RuntimeError("minimap2 not found! Install minimap2 and put it on PATH.")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"minimap2 failed:\n{e.stderr}")


def evaluate(reference, assembly_path, writer, ref_label):
    """
    Align assembly against reference, compute simple coverage/alignment metrics,
    write a CSV row.
    """
    name = os.path.basename(assembly_path)

    print(f"[*] Evaluating assembly: {assembly_path}")

    sam_text = run_minimap2(reference, assembly_path)

    total = 0
    aligned = 0
    lengths = []
    mismatches = []

    for line in sam_text.splitlines():
        if line.startswith("@"):
            continue

        total += 1
        fields = line.split("\t")

        flag = int(fields[1])
        if flag & 4:  # unmapped
            continue

        aligned += 1

        # Collect match lengths from CIGAR
        for m in re.findall(r"(\d+)M", fields[5]):
            lengths.append(int(m))

        # Collect NM:i mismatches
        for tag in fields[11:]:
            if tag.startswith("NM:i:"):
                mismatches.append(int(tag.split(":")[-1]))

    avg_len = sum(lengths) / len(lengths) if lengths else 0.0
    avg_mm = sum(mismatches) / len(mismatches) if mismatches else 0.0
    pct = (aligned / total) * 100 if total else 0.0

    writer.writerow(
        {
            "reference": ref_label,
            "assembly": name,
            "total": total,
            "aligned": aligned,
            "pct_aligned": f"{pct:.1f}",
            "avg_len": f"{avg_len:.1f}",
            "avg_mm": f"{avg_mm:.2f}",
        }
    )

    print(f"    -> {aligned}/{total} aligned "
          f"(avg_len={avg_len:.1f}, avg_mm={avg_mm:.2f})")


def append_markdown(csv_path, md_path, ref_label):
    """
    Append a Markdown summary table from the CSV file.
    """
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))

    headers = [
        "Reference",
        "Assembly",
        "Total",
        "Aligned",
        "% Aligned",
        "Avg Len",
        "Avg MM",
    ]

    table = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]

    for r in rows:
        table.append(
            f"| {r['reference']} | {r['assembly']} | {r['total']} | {r['aligned']} | {r['pct_aligned']} | {r['avg_len']} | {r['avg_mm']} |"
        )

    with open(md_path, "a") as md:
        md.write(
            textwrap.dedent(f"""
            ### Minimap2 Assembly Comparison

            **Reference genome:** `{ref_label}`

            """)
        )
        md.write("\n".join(table) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Compare multiple assembled FASTA files against a reference using minimap2"
    )

    parser.add_argument("--reference", "-r", required=True, help="Reference genome FASTA")
    parser.add_argument("--assemblies", "-a", nargs="+", required=True,
                        help="List of assembled FASTA files to evaluate")
    parser.add_argument("--output", "-o", default="assembly_eval.csv",
                        help="CSV output file")
    parser.add_argument("--md-log", "-m", default="assembly_eval.md",
                        help="Markdown summary file")

    args = parser.parse_args()

    reference = os.path.abspath(args.reference)
    if not os.path.exists(reference):
        parser.error(f"Reference not found: {reference}")
    ref_label = os.path.basename(reference)

    with open(args.output, "w", newline="") as csvfile:
        writer = csv.DictWriter(
            csvfile,
            fieldnames=[
                "reference",
                "assembly",
                "total",
                "aligned",
                "pct_aligned",
                "avg_len",
                "avg_mm",
            ],
        )
        writer.writeheader()

        for asm in args.assemblies:
            if not os.path.exists(asm):
                print(f"[skip] Assembly not found: {asm}")
                continue

            evaluate(reference, asm, writer, ref_label)

    append_markdown(args.output, args.md_log, ref_label)

    print(f"\nCSV written to: {args.output}")
    print(f"Markdown summary appended to: {args.md_log}")


if __name__ == "__main__":
    main()
