import json
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

def parse_label_from_read_id(read_id):
    read_id = str(read_id)

    if re.search(r":WGS:", read_id):
        return "WGS"

    if re.search(r":3C:", read_id):
        return "HiC"

    return "UNKNOWN"


def load_top100_insert_size_file(path):

    df = pd.read_csv(path, sep=r"\s+", compression="infer")

    required = {"read_id", "contig", "d"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Top-100 file is missing required columns: {missing}")

    out = pd.DataFrame()
    out["read_id"] = df["read_id"].astype(str)
    out["contig"] = df["contig"].astype(str)
    out["insert_size"] = pd.to_numeric(df["d"], errors="coerce")
    out["label"] = out["read_id"].apply(parse_label_from_read_id)

    out = out.dropna(subset=["insert_size", "label"])
    out = out[(out["insert_size"] > 0) & (out["label"].isin(["WGS", "HiC"]))].copy()

    return out


def em_mixture_density_log_space(params, y_min, y_max, n_grid=2500):

    y_grid = np.linspace(y_min, y_max, n_grid)
    x_grid = (10 ** y_grid) - 1

    jacobian = np.log(10) * (x_grid + 1)

    pdf_N_raw = params["pi_N"] * norm.pdf(
        x_grid,
        loc=params["mu_N"],
        scale=params["sigma_N"],
    )

    pdf_C_raw = params["pi_C"] * norm.pdf(
        x_grid,
        loc=params["mu_C"],
        scale=params["sigma_C"],
    )

    pdf_mix_raw = pdf_N_raw + pdf_C_raw

    return y_grid, pdf_mix_raw * jacobian


def plot_stacked_label_histogram(ax, insert_size, labels, bins, label_order=("WGS", "HiC")):

    x = np.asarray(insert_size, dtype=float)
    labels = np.asarray(labels)

    valid = np.isfinite(x) & (x > 0)
    x = x[valid]
    labels = labels[valid]

    y = np.log10(x + 1)

    bin_widths = np.diff(bins)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    bottom = np.zeros(len(bin_centers))

    colors = {
        "WGS": "tab:blue",
        "HiC": "tab:orange",
    }

    for lab in label_order:
        y_lab = y[labels == lab]
        counts, _ = np.histogram(y_lab, bins=bins)

        density = counts / (len(y) * bin_widths)

        ax.bar(
            bin_centers,
            density,
            width=bin_widths,
            bottom=bottom,
            align="center",
            alpha=0.45,
            color=colors.get(lab, None),
            edgecolor="none",
            label=lab,
        )

        bottom += density


def plot_top100_em_fit_by_label(
    df_top,
    params,
    title="top100",
    out_pdf=None,
    max_quantile=99.9,
):

    x = df_top["insert_size"].to_numpy(dtype=float)
    x = x[np.isfinite(x) & (x > 0)]

    if len(x) == 0:
        raise ValueError("No valid insert sizes for plotting.")

    xmax = np.percentile(x, max_quantile)
    xmax = max(xmax, params["th_N"] * 2)

    y_min = 0
    y_max = np.log10(xmax + 1)

    y_grid, mix_log = em_mixture_density_log_space(
        params=params,
        y_min=y_min,
        y_max=y_max,
    )

    y_cut = np.log10(params["th_N"] + 1)

    fig, ax = plt.subplots(figsize=(6, 4.8))
    try:
        bins_full = np.linspace(y_min, y_max, 120)

        plot_stacked_label_histogram(
            ax=ax,
            insert_size=df_top["insert_size"].to_numpy(dtype=float),
            labels=df_top["label"].to_numpy(),
            bins=bins_full,
        )

        ax.plot(
            y_grid,
            mix_log,
            color="black",
            linewidth=1.6,
            label="Fitted EM mixture distribution",
        )

        ax.axvline(
            y_cut,
            color="red",
            linestyle="--",
            linewidth=1.6,
            label="Read-selection cutoff",
        )

        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.set_xlabel(r"$\log_{10}(\mathrm{mapped\ insert\ size}+1)$")
        ax.set_ylabel("Density")
        handles, labels = ax.get_legend_handles_labels()
        handles_by_label = dict(zip(labels, handles))
        ax.legend(
            handles=[
                handles_by_label["WGS"],
                handles_by_label["HiC"],
                handles_by_label["Fitted EM mixture distribution"],
                handles_by_label["Read-selection cutoff"],
            ],
            labels=[
                "Observed shotgun insert-size distribution",
                "Observed Hi-C insert-size distribution",
                "Fitted EM mixture distribution",
                "Read-selection cutoff",
            ],
            loc="upper center",
            bbox_to_anchor=(0.5, -0.20),
            frameon=False,
            fontsize=9,
        )

        fig.tight_layout()
        if out_pdf is not None:
            fig.savefig(out_pdf, bbox_inches="tight")
    finally:
        plt.close(fig)
    

def normalize_label(x):
    x = str(x)

    if x in ["WGS", "wgs"]:
        return "WGS"

    if x in ["HiC", "HIC", "hic", "Hi-C", "3C"]:
        return "HiC"

    return "UNKNOWN"


def load_all_insert_size_label_file(path):
    df = pd.read_csv(path, sep=r"\s+", compression="infer")

    if "insert_size" not in df.columns:
        raise ValueError("All-insert-size file must contain column 'insert_size'.")

    if "label" not in df.columns:
        raise ValueError("All-insert-size file must contain column 'label'.")

    out = df[["insert_size", "label"]].copy()
    out["insert_size"] = pd.to_numeric(out["insert_size"], errors="coerce")
    out["label"] = out["label"].apply(normalize_label)

    out = out.dropna(subset=["insert_size", "label"])
    out = out[(out["insert_size"] > 0) & (out["label"].isin(["WGS", "HiC"]))].copy()

    return out


def evaluate_threshold_on_all_insert_sizes(df_all, threshold):
    df = df_all.copy()
    df["pred"] = np.where(df["insert_size"] <= threshold, "WGS", "HiC")

    tp = int(((df["label"] == "WGS") & (df["pred"] == "WGS")).sum())
    fp = int(((df["label"] == "HiC") & (df["pred"] == "WGS")).sum())
    tn = int(((df["label"] == "HiC") & (df["pred"] == "HiC")).sum())
    fn = int(((df["label"] == "WGS") & (df["pred"] == "HiC")).sum())

    precision = tp / (tp + fp) if (tp + fp) else np.nan
    recall = tp / (tp + fn) if (tp + fn) else np.nan
    specificity = tn / (tn + fp) if (tn + fp) else np.nan
    accuracy = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) else np.nan

    f1 = (
        2 * tp / (2 * tp + fp + fn)
        if (2 * tp + fp + fn) > 0
        else np.nan
    )

    em_selected = int((df["pred"] == "WGS").sum())
    em_rejected = int((df["pred"] == "HiC").sum())

    wgs_sizes = df.loc[df["label"] == "WGS", "insert_size"].to_numpy(dtype=float)
    wgs_p99 = float(np.percentile(wgs_sizes, 99))
    long_hic_cutoff = max(500, wgs_p99)

    short_hic = (df["label"] == "HiC") & (df["insert_size"] <= long_hic_cutoff)
    long_hic = (df["label"] == "HiC") & (df["insert_size"] > long_hic_cutoff)

    short_hic_retained = (
        int((short_hic & (df["pred"] == "WGS")).sum()) / int(short_hic.sum())
        if int(short_hic.sum()) > 0 else np.nan
    )

    long_hic_rejected = (
        int((long_hic & (df["pred"] == "HiC")).sum()) / int(long_hic.sum())
        if int(long_hic.sum()) > 0 else np.nan
    )

    return {
        "n_total": int(len(df)),
        "n_WGS": int((df["label"] == "WGS").sum()),
        "n_HiC": int((df["label"] == "HiC").sum()),

        "TP": tp,
        "FP": fp,
        "TN": tn,
        "FN": fn,

        "precision_WGS": precision,
        "recall_WGS": recall,
        "F1_WGS": f1,
        "accuracy": accuracy,
        "specificity_HiC": specificity,

        "EM_selected_pairs": em_selected,
        "EM_rejected_pairs": em_rejected,
        "EM_selected_fraction_of_intra_pairs": em_selected / len(df),

        "fraction_true_WGS_among_EM_selected": precision,
        "fraction_true_HiC_among_EM_selected": 1 - precision,

        "WGS_p99": wgs_p99,
        "long_HiC_cutoff": long_hic_cutoff,
        "short_insert_HiC_retained_fraction": short_hic_retained,
        "long_insert_HiC_rejected_fraction": long_hic_rejected,
    }


cases = {
    "sim3c_seed101": {
        "title": "Seed: 101",
        "top100": "sim3c_seed101.top100.insert_size.tsv.gz",
        "all": "sim3c_seed101.insert_size.label.tsv.gz",
        "em_parameters": "sim3c_seed101.em_parameters.json",
    },
}

summaries = []

for scenario, paths in cases.items():
    print(f"[INFO] Processing {scenario}")

    df_top = load_top100_insert_size_file(paths["top100"])

    with open(paths["em_parameters"], encoding="utf-8") as handle:
        params_top = json.load(handle)

    if "th_N" in params_top and "t_N" in params_top:
        if not np.isclose(float(params_top["th_N"]), float(params_top["t_N"]), rtol=1e-8, atol=1e-8):
            raise ValueError(f"Conflicting th_N and t_N in parameters for {scenario}.")
    if "th_N" not in params_top:
        params_top["th_N"] = params_top["t_N"]

    for key in ("mu_N", "sigma_N", "pi_N", "mu_C", "sigma_C", "pi_C", "th_N"):
        params_top[key] = float(params_top[key])
        if not np.isfinite(params_top[key]):
            raise ValueError(f"EM parameter {key} must be finite for {scenario}.")
    if min(params_top["sigma_N"], params_top["sigma_C"], params_top["th_N"]) <= 0:
        raise ValueError(f"Standard deviations and threshold must be positive for {scenario}.")
    if min(params_top["pi_N"], params_top["pi_C"]) < 0 or not np.isclose(
        params_top["pi_N"] + params_top["pi_C"], 1.0, rtol=1e-5, atol=1e-8,
    ):
        raise ValueError(f"Supplied mixture weights must be nonnegative and sum to one for {scenario}.")

    plot_top100_em_fit_by_label(
        df_top=df_top,
        params=params_top,
        title=paths["title"],
        out_pdf=f"{scenario}.top100_EM_fit_WGS_HiC.pdf",
    )

    df_all = load_all_insert_size_label_file(paths["all"])

    metrics = evaluate_threshold_on_all_insert_sizes(
        df_all=df_all,
        threshold=params_top["th_N"],
    )

    summary = {
        "scenario": scenario,

        "em_parameters_file": paths["em_parameters"],
        "n_top100_pairs_for_plot": len(df_top),
        "n_top100_WGS": int((df_top["label"] == "WGS").sum()),
        "n_top100_HiC": int((df_top["label"] == "HiC").sum()),
        "top100_fraction_of_all_intra_pairs": len(df_top) / len(df_all),

        "mu_N": params_top["mu_N"],
        "sigma_N": params_top["sigma_N"],
        "pi_N": params_top["pi_N"],
        "mu_C": params_top["mu_C"],
        "sigma_C": params_top["sigma_C"],
        "pi_C": params_top["pi_C"],
        "th_N": params_top["th_N"],

        **metrics,
    }

    summaries.append(summary)

summary_top100_fit_all_eval = pd.DataFrame(summaries)

summary_top100_fit_all_eval.to_csv(
    "sim3c_top100_fit_all_eval.summary.tsv",
    sep="\t",
    index=False,
)

summary_top100_fit_all_eval
