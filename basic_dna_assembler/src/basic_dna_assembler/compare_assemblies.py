import argparse
import os
import csv
import textwrap
import subprocess
import re


def evaluate_assembly(name, reference, contigs_path, writer, ref_label):
    """
    Align contigs to the reference with minimap2, parse the SAM,
    and write one CSV row of metrics including the reference label.
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
        print("ERROR: minimap2 not found! Please install minimap2.")
        return
    except subprocess.CalledProcessError as e:
        print(f"Error running minimap2 on {name}:\n{e.stderr}")
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
        f"-> {name}: {aligned}/{total} aligned (avg_len={avg_len:.1f}, avg_mm={avg_mm:.2f})"
    )


def append_markdown_table(csv_path, md_log, ref_label):
    """
    Append a Markdown summary table of the CSV results, including
    the short reference label at the top of the section.
    """
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))

    headers = [
        "Reference",
        "Assembler",
        "Total",
        "Aligned",
        "% Aligned",
        "Avg Len",
        "Avg MM",
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

        ### Assembly Comparison Results

        **Reference genome:** `{ref_label}`

        """
            )
        )
        log.write("\n".join(lines) + "\n")


def main():
    """
    Parse arguments, evaluate each assembled contigs file, and
    write both CSV metrics and a Markdown summary.
    """
    parser = argparse.ArgumentParser(
        description="Evaluate existing assembled contigs and aggregate metrics"
    )
    parser.add_argument(
        "--data-dir",
        "-d",
        default="./data",
        help="Directory containing assembled_{name}.fa files",
    )
    parser.add_argument(
        "--reference", "-r", required=True, help="Reference genome FASTA"
    )
    parser.add_argument(
        "--assemblers", "-a", nargs="+", default=["dnaasm"], help="Names of assemblers"
    )
    parser.add_argument(
        "--output",
        "-o",
        default="assembly_comparison.csv",
        help="CSV file to write metrics to",
    )
    parser.add_argument(
        "--md-log",
        "-m",
        default="comparison.log.md",
        help="Markdown log file to append summary to",
    )
    args = parser.parse_args()

    reference = os.path.abspath(args.reference)
    if not os.path.exists(reference):
        parser.error(f"Reference not found: {reference}")
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
                print(f"→ Skipping {name}: {contigs} not found")
                continue
            evaluate_assembly(name, reference, contigs, writer, ref_label)

    append_markdown_table(args.output, args.md_log, ref_label)
    print(f"\nResults written to {args.output}")
    print(f"Markdown log updated at {args.md_log}")


if __name__ == "__main__":
    main()
