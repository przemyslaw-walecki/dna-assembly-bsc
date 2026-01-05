import re
import csv
from pathlib import Path
from statistics import mean, pstdev, stdev

RUN_RE = re.compile(
    r"^\[RUN\s+(?P<run_i>\d+)/(?P<run_n>\d+)\]\s+"
    r"(?P<variant>[^|]+)\|\s*fold=(?P<fold>\d+)\s*\|\s*test=(?P<test>[^-]+)->\s*"
    r"acc_dec=(?P<acc_dec>[0-9.]+)\s+"
    r"prec=(?P<prec>[0-9.]+)\s+"
    r"rec=(?P<rec>[0-9.]+)\s+"
    r"f1=(?P<f1>[0-9.]+)\s+"
    r"acc_bub=(?P<acc_bub>[0-9.]+)\s*"
    r"\(decisions=(?P<decisions>\d+),\s*bubbles=(?P<bubbles>\d+)\)"
)

def _safe_stdev(xs):
    # klasyczne odchylenie standardowe z próby (N-1)
    return stdev(xs) if len(xs) >= 2 else 0.0

def parse_runs(log_path: Path):
    rows = []
    with log_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            m = RUN_RE.match(line)
            if not m:
                continue
            d = m.groupdict()
            rows.append({
                "variant": d["variant"].strip(),
                "fold": int(d["fold"]),
                "test": d["test"].strip(),
                "acc_dec": float(d["acc_dec"]),
                "prec": float(d["prec"]),
                "rec": float(d["rec"]),
                "f1": float(d["f1"]),
                "acc_bub": float(d["acc_bub"]),
                "decisions": int(d["decisions"]),
                "bubbles": int(d["bubbles"]),
            })
    return rows

def aggregate_all_runs(rows):
    # grupowanie: (variant) -> lista obserwacji
    by_var = {}
    for r in rows:
        by_var.setdefault(r["variant"], []).append(r)

    out = []
    for variant, xs in sorted(by_var.items()):
        acc_dec = [x["acc_dec"] for x in xs]
        prec    = [x["prec"]    for x in xs]
        rec     = [x["rec"]     for x in xs]
        f1      = [x["f1"]      for x in xs]
        acc_bub = [x["acc_bub"] for x in xs]

        out.append({
            "variant": variant,
            "n_runs_total": len(xs),
            "acc_dec_mean": mean(acc_dec),
            "acc_dec_std":  _safe_stdev(acc_dec),
            "precision_mean": mean(prec),
            "precision_std":  _safe_stdev(prec),
            "recall_mean": mean(rec),
            "recall_std":  _safe_stdev(rec),
            "f1_mean": mean(f1),
            "f1_std":  _safe_stdev(f1),
            "acc_bub_mean": mean(acc_bub),
            "acc_bub_std":  _safe_stdev(acc_bub),
        })
    return out

def aggregate_by_fold_then_std(rows):
    # grupowanie: (variant, fold) -> lista obserwacji; potem średnia w foldzie; potem std po foldach
    by_var_fold = {}
    for r in rows:
        key = (r["variant"], r["fold"])
        by_var_fold.setdefault(key, []).append(r)

    # z foldów robimy punkty (średnia z 5 runów)
    fold_points = {}
    for (variant, fold), xs in by_var_fold.items():
        fold_points.setdefault(variant, []).append({
            "fold": fold,
            "acc_dec": mean([x["acc_dec"] for x in xs]),
            "prec":    mean([x["prec"]    for x in xs]),
            "rec":     mean([x["rec"]     for x in xs]),
            "f1":      mean([x["f1"]      for x in xs]),
            "acc_bub": mean([x["acc_bub"] for x in xs]),
        })

    out = []
    for variant, pts in sorted(fold_points.items()):
        acc_dec = [p["acc_dec"] for p in pts]
        prec    = [p["prec"]    for p in pts]
        rec     = [p["rec"]     for p in pts]
        f1      = [p["f1"]      for p in pts]
        acc_bub = [p["acc_bub"] for p in pts]

        out.append({
            "variant": variant,
            "n_folds": len(pts),
            "acc_dec_mean": mean(acc_dec),
            "acc_dec_std":  _safe_stdev(acc_dec),
            "precision_mean": mean(prec),
            "precision_std":  _safe_stdev(prec),
            "recall_mean": mean(rec),
            "recall_std":  _safe_stdev(rec),
            "f1_mean": mean(f1),
            "f1_std":  _safe_stdev(f1),
            "acc_bub_mean": mean(acc_bub),
            "acc_bub_std":  _safe_stdev(acc_bub),
        })
    return out

def write_csv(path: Path, rows, fieldnames):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile", type=str, help="plik z printami (log)")
    ap.add_argument("--out_prefix", type=str, default="results_with_std",
                    help="prefiks plików wynikowych CSV")
    args = ap.parse_args()

    log_path = Path(args.logfile)
    rows = parse_runs(log_path)
    if not rows:
        raise SystemExit("Nie znaleziono żadnych linii [RUN ...] pasujących do wzorca. Sprawdź format logu.")

    allruns = aggregate_all_runs(rows)
    byfold  = aggregate_by_fold_then_std(rows)

    out1 = Path(f"{args.out_prefix}_allruns.csv")
    out2 = Path(f"{args.out_prefix}_byfold.csv")

    write_csv(out1, allruns, list(allruns[0].keys()))
    write_csv(out2, byfold,  list(byfold[0].keys()))

    print(f"Zapisano: {out1}")
    print(f"Zapisano: {out2}")

if __name__ == "__main__":
    main()