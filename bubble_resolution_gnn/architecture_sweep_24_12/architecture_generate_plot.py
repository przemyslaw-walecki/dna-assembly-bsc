import pandas as pd
import matplotlib.pyplot as plt

# pliki
base_csv = "sweep_24_12_results.csv"
agg_csv  = "architecture_sweeps_std_attached_allruns.csv"

# wczytaj i sklej params
df_base = pd.read_csv(base_csv)[["variant", "params"]]
df_agg  = pd.read_csv(agg_csv)

df = pd.merge(df_agg, df_base, on="variant", how="inner").dropna(subset=["params"])
df["params"] = df["params"].astype(int)

# wykres: tylko punkty
plt.figure(figsize=(8, 5))
plt.scatter(df["params"], df["acc_bub_mean"])

plt.xlabel("Liczba parametrów modelu")
plt.ylabel("Dokładność rozwiązywania superbąbli")
plt.title("Wpływ rozmiaru modelu na jakość rozwiązywania superbąbli")
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig("architecture_sweep_graph_acc_bub.png")