#!/usr/bin/env python3
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd

BASE = Path(__file__).resolve().parent
PROJECT = BASE.parents[1]
FONT_FAMILY = "Arial"
FONT_DIR = PROJECT / "conda_envs" / "metahit_env" / "fonts"
ARIAL_FILES = [FONT_DIR / name for name in ("arial.ttf", "arialbd.ttf", "arialbi.ttf", "ariali.ttf")]
for font_file in ARIAL_FILES:
    if not font_file.is_file():
        raise RuntimeError(f"Arial font file is missing: {font_file}")
    font_manager.fontManager.addfont(str(font_file))
font_manager.findfont(FONT_FAMILY, fallback_to_default=False)

DATASETS = [
    {
        "title": "Human",
        "before": PROJECT / "hg/6_binning/results/metahit/metahit_50_10_bins.stats",
        "after": PROJECT / "hg/7_reassembly/results/reassembly_sg_hic/reassembled_bins.stats",
    },
    {
        "title": "Pig",
        "before": PROJECT / "pig/6_binning/results/metahit/metahit_50_10_bins.stats",
        "after": PROJECT / "pig/7_reassembly/results/reassembly_sg_hic/reassembled_bins.stats",
    },
    {
        "title": "Bovine",
        "before": PROJECT / "bovine/6_binning/results/metahit/metahit_50_10_bins.stats",
        "after": PROJECT / "bovine/7_reassembly/results/reassembly_sg_hic/reassembled_bins.stats",
    },
    {
        "title": "Wastewater",
        "before": PROJECT / "ww/6_binning/results/metahit/metahit_50_10_bins.stats",
        "after": PROJECT / "ww/7_reassembly/results/reassembly_sg_hic/reassembled_bins.stats",
    },
    {
        "title": "Mats",
        "before": PROJECT / "mat/6_binning/results/metahit/metahit_50_10_bins.stats",
        "after": PROJECT / "mat/7_reassembly/results/reassembly_sg_hic/reassembled_bins.stats",
    },
]


mpl.rcParams.update(
    {
        "font.family": [FONT_FAMILY],
        "font.sans-serif": [FONT_FAMILY],
        "font.size": 12,
        "axes.titlesize": 20,
        "axes.labelsize": 12,
        "legend.fontsize": 15,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

BEFORE_COLOR = "#1f77b4"
AFTER_COLOR = "#ff7f0e"


def load_quality(path):
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t")
    required = {"completeness", "contamination"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
    return df[(df["completeness"] >= 50) & (df["contamination"] <= 10)].copy()


def matched_sorted_curves(before, after, column, ascending):
    before_sorted = before.sort_values(column, ascending=ascending).reset_index(drop=True)
    after_sorted = after.sort_values(column, ascending=ascending).reset_index(drop=True)
    n = min(len(before_sorted), len(after_sorted))
    return before_sorted.iloc[:n][column].to_numpy(), after_sorted.iloc[:n][column].to_numpy()


def plot_dataset(fig, outer_spec, dataset, show_ylabel):
    inner = outer_spec.subgridspec(2, 1, hspace=0.28)
    ax_comp = fig.add_subplot(inner[0, 0])
    ax_cont = fig.add_subplot(inner[1, 0])

    before = load_quality(dataset["before"])
    after = load_quality(dataset["after"])

    before_comp, after_comp = matched_sorted_curves(before, after, "completeness", False)
    before_cont, after_cont = matched_sorted_curves(before, after, "contamination", True)

    ax_comp.plot(before_comp, color=BEFORE_COLOR, linewidth=1.2)
    ax_comp.plot(after_comp, color=AFTER_COLOR, linewidth=1.2)
    ax_comp.set_title(dataset["title"], pad=4)
    ax_comp.set_ylim(50, 100)

    ax_cont.plot(before_cont, color=BEFORE_COLOR, linewidth=1.2)
    ax_cont.plot(after_cont, color=AFTER_COLOR, linewidth=1.2)
    ax_cont.set_ylim(0, 10)

    if show_ylabel:
        ax_comp.set_ylabel("Completion", labelpad=4)
        ax_cont.set_ylabel("Contamination", labelpad=4)

    for ax in (ax_comp, ax_cont):
        ax.spines["top"].set_color("#999999")
        ax.spines["right"].set_color("#999999")
        ax.spines["left"].set_color("#999999")
        ax.spines["bottom"].set_color("#999999")
        ax.tick_params(length=3, width=0.8)


def main():
    fig = plt.figure(figsize=(11.5, 7.4), facecolor="white")
    gs = fig.add_gridspec(
        2,
        6,
        left=0.08,
        right=0.98,
        top=0.91,
        bottom=0.19,
        wspace=0.40,
        hspace=0.48,
    )

    positions = [
        gs[0, 0:2],
        gs[0, 2:4],
        gs[0, 4:6],
        gs[1, 1:3],
        gs[1, 3:5],
    ]

    for i, (dataset, pos) in enumerate(zip(DATASETS, positions)):
        plot_dataset(fig, pos, dataset, show_ylabel=i in (0, 3))

    fig.text(0.52, 0.13, "Bin Ranking", ha="center", va="center", fontsize=13)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=BEFORE_COLOR, alpha=0.8, label="Before reassembly"),
        plt.Rectangle((0, 0), 1, 1, color=AFTER_COLOR, alpha=0.8, label="After reassembly"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.04))

    for ext in ("png", "pdf"):
        fig.savefig(BASE / f"bin_quality.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {BASE / 'bin_quality.png'}")
    print(f"Wrote {BASE / 'bin_quality.pdf'}")


if __name__ == "__main__":
    main()
