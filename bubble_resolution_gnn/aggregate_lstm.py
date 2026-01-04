import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

IN_CSV =  "./lstm_vs_mean/lstm_vs_mean_cov20.csv"

# --- Wczytanie
df = pd.read_csv(IN_CSV)

# Upewnij się, że kluczowe kolumny istnieją
required_cols = ["use_lstm","fold","test_genome","acc_bub","acc_dec","brier","ece"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise KeyError(f"Brak kolumn w pliku wejściowym: {missing}")

# --- Agregacja globalna (po wszystkich foldach i genomach)
agg_global = (
    df.groupby("use_lstm")
      .agg(acc_bub_mean=("acc_bub","mean"),
           acc_bub_std=("acc_bub","std"),
           acc_dec_mean=("acc_dec","mean"),
           acc_dec_std=("acc_dec","std"),
           brier_mean=("brier","mean"),
           brier_std=("brier","std"),
           ece_mean=("ece","mean"),
           ece_std=("ece","std"),
           runs=("acc_bub","count"))
      .reset_index()
)

# Mapowanie na etykiety PL
label_map = {False: "Bez LSTM", True: "Z LSTM"}
agg_global["wariant"] = agg_global["use_lstm"].map(label_map)

# Zaokrąglenia do 2 miejsc
for c in ["acc_bub_mean","acc_bub_std","acc_dec_mean","acc_dec_std","brier_mean","brier_std","ece_mean","ece_std"]:
    agg_global[c] = agg_global[c].round(2)

# --- Agregacja per genom (dla tabeli w pracy)
agg_genome = (
    df.groupby(["test_genome","use_lstm"])
      .agg(acc_bub_mean=("acc_bub","mean"),
           acc_bub_std=("acc_bub","std"),
           acc_dec_mean=("acc_dec","mean"),
           acc_dec_std=("acc_dec","std"),
           brier_mean=("brier","mean"),
           brier_std=("brier","std"),
           ece_mean=("ece","mean"),
           ece_std=("ece","std"),
           runs=("acc_bub","count"))
      .reset_index()
)
agg_genome["wariant"] = agg_genome["use_lstm"].map(label_map)
for c in ["acc_bub_mean","acc_bub_std","acc_dec_mean","acc_dec_std","brier_mean","brier_std","ece_mean","ece_std"]:
    agg_genome[c] = agg_genome[c].round(2)

# --- Zapis CSV
agg_global[["wariant","runs","acc_bub_mean","acc_bub_std","acc_dec_mean","acc_dec_std","brier_mean","brier_std","ece_mean","ece_std"]] \
    .to_csv("podsumowanie_lstm_cov20.csv", index=False, encoding="utf-8")

agg_genome[["test_genome","wariant","runs","acc_bub_mean","acc_bub_std","acc_dec_mean","acc_dec_std","brier_mean","brier_std","ece_mean","ece_std"]] \
    .to_csv("podsumowanie_lstm_cov20_per_genom.csv", index=False, encoding="utf-8")

# --- Wykres 1: Dokładność (acc_bub, acc_dec) z błędami
x = np.arange(len(agg_global))
width = 0.35

plt.figure(figsize=(7,4))
plt.bar(x - width/2, agg_global["acc_bub_mean"], width, yerr=agg_global["acc_bub_std"], capsize=4, label="Dokładność bąbli")
plt.bar(x + width/2, agg_global["acc_dec_mean"], width, yerr=agg_global["acc_dec_std"], capsize=4, label="Dokładność decyzji")
plt.xticks(x, agg_global["wariant"])
plt.ylabel("Dokładność")
plt.title("Porównanie LSTM vs bez LSTM (pokrycie 20×)")
plt.legend()
plt.tight_layout()
plt.savefig("porownanie_lstm_cov20_dokladnosc.png", dpi=200)
plt.close()

# --- Wykres 2: Kalibracja (Brier, ECE) z błędami
plt.figure(figsize=(7,4))
plt.bar(x - width/2, agg_global["brier_mean"], width, yerr=agg_global["brier_std"], capsize=4, label="Błąd Brier’a")
plt.bar(x + width/2, agg_global["ece_mean"],   width, yerr=agg_global["ece_std"],   capsize=4, label="ECE")
plt.xticks(x, agg_global["wariant"])
plt.ylabel("Wartość")
plt.title("Kalibracja: LSTM vs bez LSTM")
plt.legend()
plt.tight_layout()
plt.savefig("porownanie_lstm_cov20_kalibracja.png", dpi=200)
plt.close()

print("Zapisano: podsumowanie_lstm_cov20.csv, podsumowanie_lstm_cov20_per_genom.csv,")
print("          porownanie_lstm_cov20_dokladnosc.png, porownanie_lstm_cov20_kalibracja.png")
