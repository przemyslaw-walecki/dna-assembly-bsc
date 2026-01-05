import matplotlib.pyplot as plt

# dane z tabeli
coverage = [10, 20, 30]
acc_dec = [0.7114624505928854, 0.5884805837894188, 0.5078763127187864]

plt.figure(figsize=(7, 4))
plt.plot(coverage, acc_dec, marker="o", linestyle="-")

plt.xlabel("Pokrycie [x]")
plt.ylabel("Dokładność")
plt.title("Wpływ poziomu pokrycia na dokładność rozwiązywania superbąbli")

plt.xticks(coverage, [f"{c}x" for c in coverage])
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()

plt.savefig("bub_dec_vs_coverage.png")