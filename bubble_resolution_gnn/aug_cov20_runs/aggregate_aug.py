import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

IN_CSV = "aug_vs_plain_cov20.csv"       # ścieżka do wejściowego pliku
OUT_DIR = Path("aug_vs_plain_cov20_out")
OUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(IN_CSV)

# --- agregacja globalna (po wszystkich genomach) ---

global_agg = (
    df.groupby("augmented_train")
      .agg(
          acc_bub_mean=("acc_bub", "mean"),
          acc_dec_mean=("acc_dec", "mean"),
          mrr_mean=("mrr", "mean"),
          brier_mean=("brier", "mean"),
          ece_mean=("ece", "mean"),
      )
      .reset_index()
)

# --- agregacja per genom ---

per_genome = (
    df.groupby(["test_genome", "augmented_train"])
      .agg(
          acc_bub_mean=("acc_bub", "mean"),
          acc_dec_mean=("acc_dec", "mean"),
          mrr_mean=("mrr", "mean"),
          brier_mean=("brier", "mean"),
          ece_mean=("ece", "mean"),
      )
      .reset_index()
)

# zaokrąglenie do 2 miejsc na potrzeby pracy
for c in ["acc_bub_mean", "acc_dec_mean", "mrr_mean", "brier_mean", "ece_mean"]:
    global_agg[c] = global_agg[c].round(2)
    per_genome[c] = per_genome[c].round(2)

# zapis pomocniczych CSV
global_agg.to_csv(OUT_DIR / "summary_aug_global.csv", index=False)
per_genome.to_csv(OUT_DIR / "summary_aug_per_genome.csv", index=False)

# --- wykresy: globalne porównanie augmentacja vs bez augmentacji ---


# uporządkowanie: False (bez) -> True (z)
g = global_agg.sort_values("augmented_train").reset_index(drop=True)
labels = ["Bez augmentacji", "Z augmentacją"]
x = np.arange(len(labels))
width = 0.35

# Dokładność: acc_bub vs acc_dec
fig, ax = plt.subplots(figsize=(6, 4))
bars_bub = ax.bar(x - width/2, g["acc_bub_mean"], width, label="Dokładność bąbli")
bars_dec = ax.bar(x + width/2, g["acc_dec_mean"], width, label="Dokładność decyzji")

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("Dokładność")
ax.set_title("Porównanie trenowania z augmentacją i bez (pokrycie 20×)")
ax.set_ylim(0.0, 1.0)

# etykiety nad słupkami
ax.bar_label(bars_bub, fmt="%.2f", padding=3)
ax.bar_label(bars_dec, fmt="%.2f", padding=3)

# legenda poza obszarem rysunku (po prawej)
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.)

fig.tight_layout()
fig.savefig(OUT_DIR / "porownanie_aug_cov20_dokladnosc.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# Kalibracja: Brier vs ECE
fig, ax = plt.subplots(figsize=(6, 4))
bars_brier = ax.bar(x - width/2, g["brier_mean"], width, label="Błąd Briera")
bars_ece   = ax.bar(x + width/2, g["ece_mean"],   width, label="ECE")

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("Wartość miary")
ax.set_title("Porównanie kalibracji z augmentacją i bez (pokrycie 20×)")

# etykiety nad słupkami
ax.bar_label(bars_brier, fmt="%.2f", padding=3)
ax.bar_label(bars_ece,   fmt="%.2f", padding=3)

# legenda poza obszarem rysunku
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.)

fig.tight_layout()
fig.savefig(OUT_DIR / "porownanie_aug_cov20_kalibracja.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print("Zaktualizowane wykresy zapisane w:", OUT_DIR)