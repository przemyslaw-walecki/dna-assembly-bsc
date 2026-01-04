# %% Analiza odporności modelu GNN na poziom pokrycia

import pandas as pd
import matplotlib.pyplot as plt

# Wczytanie wyników testów dla różnych poziomów pokrycia
df_cov = pd.read_csv("./cov_sweep_runs/cov_all_results.csv")

# Obliczenie średnich i odchyleń standardowych dokładności
cov_summary = (
    df_cov.groupby("coverage")
          .agg(acc_bub_mean=("acc_bub","mean"),
               acc_bub_std=("acc_bub","std"),
               acc_dec_mean=("acc_dec","mean"),
               acc_dec_std=("acc_dec","std"))
          .reset_index()
          .round(2)
)

# Zapis podsumowania do pliku CSV
cov_summary.to_csv("podsumowanie_pokrycia.csv", index=False)

# Wykres — dokładność modelu względem poziomu pokrycia
plt.figure(figsize=(6,4))
plt.errorbar(cov_summary.coverage, cov_summary.acc_bub_mean,
             yerr=cov_summary.acc_bub_std, label="Dokładność rozwiązywania bąbli", marker='o')
plt.errorbar(cov_summary.coverage, cov_summary.acc_dec_mean,
             yerr=cov_summary.acc_dec_std, label="Dokładność decyzji", marker='s')

plt.xlabel("Poziom pokrycia (coverage)")
plt.ylabel("Dokładność")
plt.title("Odporność modelu GNN na różne poziomy pokrycia")
plt.legend()
plt.tight_layout()
plt.savefig("odpornosc_na_pokrycie.png", dpi=200)
plt.close()

print("Zapisano: podsumowanie_pokrycia.csv, odpornosc_na_pokrycie.png")
