# %% Wizualizacja wyników – podział na dwa wykresy + dokładność względem liczby parametrów

import pandas as pd
import matplotlib.pyplot as plt

# Wczytanie zsumowanej tabeli wyników
df = pd.read_csv("summary_architectures.csv")  # kolumny: variant, params, acc_bub_mean, acc_dec_mean, ...

# Uporządkowanie nazw wariantów (usunięcie "_mean")
df["variant_clean"] = df["variant"].str.replace("_mean", "", regex=False)

# ---------- Podział na dwa wykresy względem mediany liczby parametrów ----------
median_params = df["params"].median()
small = df[df["params"] <= median_params].sort_values("acc_bub_mean", ascending=False).copy()
large = df[df["params"] >  median_params].sort_values("acc_bub_mean", ascending=False).copy()

def plot_bar(sub, title, outpng):
    ax = sub.set_index("variant_clean")[["acc_bub_mean","acc_dec_mean"]].plot(kind="bar", rot=45, figsize=(11,4))
    ax.set_ylabel("Dokładność")
    ax.set_title(title)
    ax.legend(["Dokładność rozwiązywania bąbli", "Dokładność decyzji"])
    ax.figure.tight_layout()
    ax.figure.savefig(outpng, dpi=200)
    plt.close(ax.figure)

plot_bar(small,
         f"Porównanie architektur (≤ mediana liczby parametrów = {int(median_params)})",
         "porownanie_architektur_male.png")

plot_bar(large,
         f"Porównanie architektur (> mediana liczby parametrów = {int(median_params)})",
         "porownanie_architektur_duze.png")

# ---------- Wykres zależności dokładności od liczby parametrów ----------
by_params = (df.groupby("params", as_index=False)
               .agg(acc_bub_mean=("acc_bub_mean","mean"),
                    acc_dec_mean=("acc_dec_mean","mean"))
               .sort_values("params"))

plt.figure(figsize=(7,4))
plt.plot(by_params["params"], by_params["acc_bub_mean"], marker="o", label="Dokładność bąbli")
plt.plot(by_params["params"], by_params["acc_dec_mean"], marker="s", label="Dokładność decyzji")
plt.xlabel("Liczba parametrów modelu")
plt.ylabel("Dokładność")
plt.title("Zależność dokładności od liczby parametrów modelu")
plt.legend()
plt.tight_layout()
plt.savefig("dokladnosc_vs_parametry.png", dpi=200)
plt.close()

print("Zapisano: porownanie_architektur_male.png, porownanie_architektur_duze.png, dokladnosc_vs_parametry.png")
