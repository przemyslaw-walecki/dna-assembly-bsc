# %% [SELECT BEST MODEL FROM SWEEP CSV]
import pandas as pd

# Load your sweep file (edit path if needed)
df = pd.read_csv("sweep_results.csv")  # or cov_all_results.csv

# Sort by priority: bubble acc ↓, decision acc ↓, brier ↑ (smaller better)
df_sorted = (
    df.sort_values(by=["acc_bub", "acc_dec", "brier"], ascending=[False, False, True])
    .reset_index(drop=True)
)

best_row = df_sorted.iloc[0]
best_variant = best_row["variant"]
print("=== Best model setup ===")
print(best_row[["variant", "acc_bub", "acc_dec", "brier", "params"]])

# Save this variant name for reuse
with open("best_variant.txt", "w") as f:
    f.write(best_variant)
print(f"Saved best variant to best_variant.txt")
