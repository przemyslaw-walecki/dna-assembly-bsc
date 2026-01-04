#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

RESULTS_CSV = Path("cov_sweep_runs") / "cov_all_results.csv"
OUT_CSV     = Path("cov_sweep_runs") / "best_models_by_genome_cov20.csv"
CKPT_DIR    = Path("cov_sweep_runs")  # tam gdzie zapisywałeś .pth

COVERAGE_TARGET = 20  # szukamy tylko modeli dla pokrycia 20×

def build_ckpt_path(row) -> str:
    """
    Odtwarza nazwę pliku checkpointu na podstawie coverage, fold, test_genome i variant.
    Zakłada schemat: cov{cov}_fold{fold}_{genome}_{variant}.pth
    """
    try:
        cov  = int(row["coverage"])
        fold = int(row["fold"])
        genome  = str(row["test_genome"])
        variant = str(row["variant"])
    except Exception:
        return ""
    tag = f"cov{cov}_fold{fold}_{genome}_{variant}"
    return str((CKPT_DIR / f"{tag}.pth").as_posix())

def main():
    if not RESULTS_CSV.is_file():
        raise FileNotFoundError(f"Nie znaleziono pliku z wynikami: {RESULTS_CSV}")

    df = pd.read_csv(RESULTS_CSV)

    required = {"test_genome", "variant", "acc_bub", "acc_dec", "coverage", "fold"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Brak wymaganych kolumn w {RESULTS_CSV}: {missing}")

    # upewniamy się, że coverage jest liczbą i filtrujemy tylko 20×
    df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce")
    df20 = df[df["coverage"] == COVERAGE_TARGET].copy()
    if df20.empty:
        raise ValueError(f"Brak wyników dla coverage == {COVERAGE_TARGET} w {RESULTS_CSV}")

    # sort: w obrębie każdego genomu najlepsze (acc_bub, acc_dec) na górze
    df20_sorted = df20.sort_values(
        by=["test_genome", "acc_bub", "acc_dec"],
        ascending=[True, False, False]
    )

    # pierwszy wiersz per genom = najlepszy wariant przy pokryciu 20×
    best = df20_sorted.groupby("test_genome", as_index=False).first()

    # odtworzenie ścieżki do checkpointa
    best["ckpt_path"] = best.apply(build_ckpt_path, axis=1)

    # zapis całości
    best.to_csv(OUT_CSV, index=False)
    print(f"Zapisano zestawienie najlepszych modeli (coverage={COVERAGE_TARGET}) do: {OUT_CSV}")

    # podgląd
    cols = ["test_genome", "coverage", "fold", "variant", "acc_bub", "acc_dec", "ckpt_path"]
    if "params" in best.columns:
        cols.append("params")

    print("\nNajlepszy model dla każdego genomu przy pokryciu 20× (gdy był testowy):\n")
    print(best[cols].to_string(index=False))


if __name__ == "__main__":
    main()
