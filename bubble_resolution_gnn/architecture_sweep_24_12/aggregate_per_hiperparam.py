# aggregate_grid_effects.py
# Usage:
#   python aggregate_grid_effects.py sweep_24_12_results.csv
#
# Output:
#   - prints 3 pretty tables (for ch / g / e) to stdout
#   - saves 3 CSV files next to input:
#       agg_by_ch.csv, agg_by_g.csv, agg_by_e.csv
#
# Aggregation:
#   For each hyperparameter value (e.g. ch=32), averages metrics across ALL
#   combinations of the other two hyperparameters (and across whatever is already
#   aggregated in your CSV). No std.

from __future__ import annotations
import re
import sys
from pathlib import Path

import pandas as pd

VAR_RE = re.compile(
    r"^(?P<ch>(?:emb|ch)\d+)_g(?P<g>\d+)_e(?P<e>\d+)(?:_mean)?$",
    re.IGNORECASE
)

# Polish column names (no abbreviations in the output tables)
POLISH_COLS = {
    "params": "Liczba parametrów",
    "acc_dec": "Dokładność decyzji",
    "precision_dec": "Precyzja decyzji",
    "recall_dec": "Czułość decyzji",
    "f1_dec": "Miara F1 decyzji",
    "acc_bub": "Dokładność superbąbli",
}

# Order of metrics in tables
METRIC_COLS = ["acc_dec", "precision_dec", "recall_dec", "f1_dec", "acc_bub", "params"]


def parse_variant(variant: str) -> tuple[int, int, int]:
    v = variant.strip()
    m = VAR_RE.match(v)
    if not m:
        raise ValueError(f"Nie rozpoznano formatu wariantu: {variant!r}")
    ch_raw = m.group("ch").lower()
    # emb8 -> ch8
    ch = int(ch_raw.replace("emb", "").replace("ch", ""))
    g = int(m.group("g"))
    e = int(m.group("e"))
    return ch, g, e


def add_hparams(df: pd.DataFrame) -> pd.DataFrame:
    ch_vals, g_vals, e_vals = [], [], []
    for v in df["variant"].astype(str).tolist():
        ch, g, e = parse_variant(v)
        ch_vals.append(ch)
        g_vals.append(g)
        e_vals.append(e)
    out = df.copy()
    out["ch"] = ch_vals
    out["g"] = g_vals
    out["e"] = e_vals
    return out


def agg_main_effect(df: pd.DataFrame, key: str) -> pd.DataFrame:
    # Mean effect by one hyperparameter, averaged over the remaining two.
    # Also keep mean params (handy as a "typical" size for that hyperparam value).
    g = df.groupby(key, as_index=False)[METRIC_COLS].mean(numeric_only=True)
    g = g.sort_values(key).reset_index(drop=True)
    return g


def format_table(df: pd.DataFrame, key: str) -> pd.DataFrame:
    # Rename columns to Polish and format numeric values nicely.
    out = df.copy()
    out = out.rename(columns={key: {"ch": "Kanały CNN (ch)", "g": "Wymiar GCN (g)", "e": "Wymiar klasyfikatora (e)"}[key]})
    out = out.rename(columns=POLISH_COLS)

    # round metrics; keep params as int-ish
    for c in out.columns:
        if c == "Liczba parametrów":
            out[c] = out[c].round(0).astype(int)
        elif c != out.columns[0]:
            out[c] = out[c].round(6)
    return out


def main():
    if len(sys.argv) < 2:
        print("Podaj ścieżkę do sweep_24_12_results.csv", file=sys.stderr)
        sys.exit(1)

    in_path = Path(sys.argv[1])
    df = pd.read_csv(in_path)

    # Basic validation
    required = {"variant", "params", "acc_dec", "precision_dec", "recall_dec", "f1_dec", "acc_bub"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Brak wymaganych kolumn w CSV: {sorted(missing)}")

    df = add_hparams(df)

    # Aggregate main effects
    by_ch = agg_main_effect(df, "ch")
    by_g  = agg_main_effect(df, "g")
    by_e  = agg_main_effect(df, "e")

    # Format for easy copy/paste (markdown-friendly)
    t_ch = format_table(by_ch, "ch")
    t_g  = format_table(by_g, "g")
    t_e  = format_table(by_e, "e")

    # Print as markdown tables (easy to paste into notes / convert to LaTeX later)
    pd.set_option("display.max_columns", 1000)
    pd.set_option("display.width", 200)

    print("\n=== Wpływ hiperparametru: kanały CNN (ch) — uśrednione po pozostałych ===\n")
    print(t_ch.to_markdown(index=False))

    print("\n\n=== Wpływ hiperparametru: wymiar GCN (g) — uśrednione po pozostałych ===\n")
    print(t_g.to_markdown(index=False))

    print("\n\n=== Wpływ hiperparametru: wymiar klasyfikatora (e) — uśrednione po pozostałych ===\n")
    print(t_e.to_markdown(index=False))

    # Save CSVs (polish headers too, so you can import directly)
    out_dir = in_path.parent
    t_ch.to_csv(out_dir / "agg_by_ch.csv", index=False)
    t_g.to_csv(out_dir / "agg_by_g.csv", index=False)
    t_e.to_csv(out_dir / "agg_by_e.csv", index=False)

    print("\n\nZapisano pliki:")
    print(f" - {out_dir / 'agg_by_ch.csv'}")
    print(f" - {out_dir / 'agg_by_g.csv'}")
    print(f" - {out_dir / 'agg_by_e.csv'}")


if __name__ == "__main__":
    main()
