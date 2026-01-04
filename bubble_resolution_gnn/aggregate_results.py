import pandas as pd
from matplotlib import pyplot as plt

df = pd.read_csv("./sweep_from_paths/sweep_results.csv")
agg = (df.groupby("variant")
       .agg(acc_bub_mean=("acc_bub","mean"),
            acc_dec_mean=("acc_dec","mean"),
            mrr_mean=("mrr","mean"),
            brier_mean=("brier","mean"),
            params=("params","first"))
       .reset_index()
       .sort_values("acc_bub_mean", ascending=False)
       .round(2))
agg.to_csv("summary_architectures.csv", index=False)
