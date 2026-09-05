import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm


dataset = "ww"
plot_title = "Wastewater"
insert_size_file = "ww.hic.top100.insert_size.tsv.gz"
em_parameters_file = "ww.em_parameters.json"
out_pdf = "ww.top100_EM_fit.pdf"

max_quantile = 99.9


arrays = []
for chunk in pd.read_csv(
    insert_size_file, sep="\t", usecols=["d"], compression="infer", chunksize=2_000_000,
):
    values = pd.to_numeric(chunk["d"], errors="coerce").to_numpy(dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if len(values):
        arrays.append(values)

x = np.concatenate(arrays)
del arrays

with open(em_parameters_file, encoding="utf-8") as handle:
    params = json.load(handle)

mu_N = float(params["mu_N"])
sigma_N = float(params["sigma_N"])
pi_N = float(params["pi_N"])
mu_C = float(params["mu_C"])
sigma_C = float(params["sigma_C"])
pi_C = float(params["pi_C"])
threshold = float(params["th_N"] if "th_N" in params else params["t_N"])


print(f"Plotting {len(x):,} insert sizes")


xmax = max(float(np.percentile(x, max_quantile)), threshold * 2)
y = np.log10(x + 1)
y_cut = np.log10(threshold + 1)
ymax = np.log10(xmax + 1)
edges = np.linspace(0, ymax, 121)

counts, _ = np.histogram(y, bins=edges)
density = counts / (len(x) * np.diff(edges))

y_grid = np.linspace(0, ymax, 2500)
x_grid = 10 ** y_grid - 1
jacobian = np.log(10) * (x_grid + 1)
mixture = (
    pi_N * norm.pdf(x_grid, loc=mu_N, scale=sigma_N)
    + pi_C * norm.pdf(x_grid, loc=mu_C, scale=sigma_C)
) * jacobian

fig, ax = plt.subplots(figsize=(6, 4.8))
try:
    observed = ax.bar(
        edges[:-1], density, width=np.diff(edges), align="edge",
        alpha=0.45, color="tab:blue", linewidth=0,
        label="Observed insert-size distribution",
    )
    fitted, = ax.plot(
        y_grid, mixture, color="black", linewidth=1.2,
        label="Fitted EM mixture distribution",
    )
    cutoff = ax.axvline(
        y_cut, color="red", linestyle="--", linewidth=1.2,
        label="Read-selection cutoff",
    )
    ax.set_xlim(0, ymax)
    ax.set_ylim(bottom=0)
    ax.set_title(plot_title, fontsize=14, fontweight="bold")
    ax.set_xlabel(r"$\log_{10}(\mathrm{mapped\ insert\ size}+1)$")
    ax.set_ylabel("Density")
    ax.legend(
        handles=[observed, fitted, cutoff], loc="upper center",
        bbox_to_anchor=(0.5, -0.20), frameon=False, fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
finally:
    plt.close(fig)

print(f"Saved: {out_pdf}")
