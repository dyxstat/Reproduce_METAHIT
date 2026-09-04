#!/usr/bin/env python3
"""Calculate and plot supplementary taxonomic breadth from current results."""

from __future__ import annotations

import os
import re
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


FONT_FAMILY = "Arial"
FONT_DIR = ROOT / "conda_envs" / "metahit_env" / "fonts"
ARIAL_FILES = [FONT_DIR / name for name in ("arial.ttf", "arialbd.ttf", "arialbi.ttf", "ariali.ttf")]
for font_file in ARIAL_FILES:
    if not font_file.is_file():
        raise RuntimeError(f"Arial font file is missing: {font_file}")
    font_manager.fontManager.addfont(str(font_file))
font_manager.findfont(FONT_FAMILY, fallback_to_default=False)

REFERENCE = "METAHICT"
LEVELS = ["Species", "Genus", "Family", "Order"]
COMPARISONS = ["bin3C", "MetaCC", "ImputeCC"]

COLORS = {
    "Both": "#7f8c8d",
    f"{REFERENCE} only": "#3498db",
    "bin3C only": "#9b59b6",
    "MetaCC only": "#f39c12",
    "ImputeCC only": "#e74c3c",
}

DATASET_INPUTS = {
    "Sheep Gut": {
        "dataset": "sheep",
        "metahict_annotation": ROOT / "sheep/9_annotation/results_binners/annotation_metahit/classify",
    },
    "Cow Rumen": {
        "dataset": "cow",
        "metahict_annotation": ROOT / "cow/9_annotation/results/classify",
    },
    "Bovine Skin": {
        "dataset": "bovine",
        "metahict_annotation": ROOT / "bovine/9_annotation/results_binners/annotation_metahit/classify",
    },
    "Hydrothermal Mats": {
        "dataset": "mat",
        "metahict_annotation": ROOT / "mat/9_annotation/results_binners/annotation_metahit/classify",
    },
}

DATA: dict[str, dict[str, dict[str, list[int]]]] = {}


def load_annotations(annotation_dir: Path) -> dict[str, str]:
    annotations: dict[str, str] = {}
    for marker in ("bac120", "ar53"):
        summary = annotation_dir / f"gtdbtk.{marker}.summary.tsv"
        if not summary.exists():
            continue
        frame = pd.read_csv(summary, sep="\t")
        for bin_id, classification in zip(
            frame["user_genome"].astype(str), frame["classification"].astype(str)
        ):
            normalized = re.sub(r"\\.(?:orig|strict|permissive)$", "", bin_id)
            normalized = re.sub(r"^bin(\\d+)$", r"bin.\\1", normalized)
            annotations[normalized] = classification
    if not annotations:
        raise FileNotFoundError(f"No GTDB-Tk summary files found in {annotation_dir}")
    return annotations


def load_medium_quality_bins(stats_file: Path) -> set[str]:
    if not stats_file.exists():
        raise FileNotFoundError(f"Missing bin-quality table: {stats_file}")
    frame = pd.read_csv(stats_file, sep="\t")
    selected = frame[(frame["completeness"] >= 50) & (frame["contamination"] < 10)]
    return set(selected.iloc[:, 0].astype(str))


def taxonomic_sets(annotation_dir: Path, stats_file: Path) -> dict[str, set[str]]:
    annotations = load_annotations(annotation_dir)
    selected_bins = load_medium_quality_bins(stats_file)
    taxa = {level: set() for level in LEVELS}
    prefixes = {"Species": "s__", "Genus": "g__", "Family": "f__", "Order": "o__"}
    for bin_id in selected_bins:
        classification = annotations.get(bin_id)
        if not classification:
            continue
        fields = str(classification).split(";")
        for level, prefix in prefixes.items():
            value = next((field[len(prefix):] for field in fields if field.startswith(prefix)), "")
            if value:
                taxa[level].add(value)
    return taxa


def calculate_current_data() -> dict[str, dict[str, dict[str, list[int]]]]:
    output = {}
    for dataset_label, config in DATASET_INPUTS.items():
        dataset = config["dataset"]
        result_root = ROOT / dataset / "6_binning/results/metahit"
        annotation_root = ROOT / dataset / "9_annotation/results_binners"
        annotation_dirs = {
            REFERENCE: config["metahict_annotation"],
            "bin3C": annotation_root / "annotation_bin3c/classify",
            "MetaCC": annotation_root / "annotation_metacc/classify",
            "ImputeCC": annotation_root / "annotation_imputecc/classify",
        }
        stats_files = {
            REFERENCE: result_root / "metahit_50_10_bins.stats",
            "bin3C": result_root / "work_files/binsC.stats",
            "MetaCC": result_root / "work_files/binsA.stats",
            "ImputeCC": result_root / "work_files/binsB.stats",
        }
        taxa = {
            tool: taxonomic_sets(annotation_dirs[tool], stats_files[tool])
            for tool in (REFERENCE, *COMPARISONS)
        }
        output[dataset_label] = {}
        for comparison in COMPARISONS:
            output[dataset_label][comparison] = {
                "Both": [len(taxa[REFERENCE][level] & taxa[comparison][level]) for level in LEVELS],
                f"{REFERENCE} only": [
                    len(taxa[REFERENCE][level] - taxa[comparison][level]) for level in LEVELS
                ],
                f"{comparison} only": [
                    len(taxa[comparison][level] - taxa[REFERENCE][level]) for level in LEVELS
                ],
            }
    return output


DATA = calculate_current_data()

def configure_style() -> None:
    plt.rcParams.update({
        "font.family": [FONT_FAMILY],
        "font.sans-serif": [FONT_FAMILY],
        "font.size": 30,
        "axes.labelsize": 30,
        "axes.titlesize": 35,
        "xtick.labelsize": 30,
        "ytick.labelsize": 30,
        "legend.fontsize": 24,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def add_labels(ax, x_pos, both_values, ref_only_values, comp_only_values, split_width):
    for idx in range(len(LEVELS)):
        both = both_values[idx]
        ref_only = ref_only_values[idx]
        comp_only = comp_only_values[idx]

        if both > 0:
            ax.text(x_pos[idx], both / 2, str(both), ha="center", va="center", color="white", fontsize=25)
        if comp_only > 0:
            ax.text(
                x_pos[idx] - split_width / 2,
                both + comp_only + 0.5,
                str(comp_only),
                ha="center",
                va="bottom",
                color="black",
                fontsize=25,
            )
        if ref_only > 0:
            ax.text(
                x_pos[idx] + split_width / 2,
                both + ref_only + 0.5,
                str(ref_only),
                ha="center",
                va="bottom",
                color="black",
                fontsize=25,
            )


def draw_dataset_row(axes, dataset: str) -> None:
    bar_width = 0.6
    split_width = bar_width / 2
    x_pos = np.arange(len(LEVELS))

    for col_idx, comp_tool in enumerate(COMPARISONS):
        ax = axes[col_idx]
        comp_only_label = f"{comp_tool} only"
        values = DATA[dataset][comp_tool]
        both_values = np.array(values["Both"])
        ref_only_values = np.array(values[f"{REFERENCE} only"])
        comp_only_values = np.array(values[comp_only_label])

        ax.bar(x_pos, both_values, bar_width, color=COLORS["Both"], label="Both", alpha=0.8, edgecolor="black")
        ax.bar(
            x_pos - split_width / 2,
            comp_only_values,
            split_width,
            bottom=both_values,
            color=COLORS[comp_only_label],
            label=comp_only_label,
            alpha=0.8,
            edgecolor="black",
        )
        ax.bar(
            x_pos + split_width / 2,
            ref_only_values,
            split_width,
            bottom=both_values,
            color=COLORS[f"{REFERENCE} only"],
            label=f"{REFERENCE} only",
            alpha=0.8,
            edgecolor="black",
        )

        add_labels(ax, x_pos, both_values, ref_only_values, comp_only_values, split_width)

        ax.set_ylabel("Number of taxa", fontsize=30)
        ax.set_title(f"{REFERENCE} vs {comp_tool}", fontsize=35)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(LEVELS, fontsize=30)
        ax.grid(True, alpha=0.3, axis="y")
        ax.set_axisbelow(True)

        max_height = int(np.max(both_values + np.maximum(ref_only_values, comp_only_values)))
        ax.set_ylim(0, max_height * 1.15 if max_height > 0 else 10)


def plot() -> None:
    configure_style()
    datasets = list(DATA)
    fig, axes = plt.subplots(len(datasets), len(COMPARISONS), figsize=(24, 32), sharex=False, sharey=False)
    axes = np.atleast_2d(axes)

    for row_idx, dataset in enumerate(datasets):
        draw_dataset_row(axes[row_idx], dataset)
        axes[row_idx, 0].text(
            -0.23,
            1.04,
            chr(ord("A") + row_idx),
            transform=axes[row_idx, 0].transAxes,
            fontsize=45,
            fontweight="normal",
            ha="right",
            va="top",
        )

    legend_handles = [
        Patch(facecolor=COLORS["Both"], label="Both"),
        Patch(facecolor=COLORS[f"{REFERENCE} only"], label=f"{REFERENCE} only"),
        Patch(facecolor=COLORS["bin3C only"], label="bin3C only"),
        Patch(facecolor=COLORS["MetaCC only"], label="MetaCC only"),
        Patch(facecolor=COLORS["ImputeCC only"], label="ImputeCC only"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.5, 0.01),
        columnspacing=1.8,
        handlelength=2.0,
    )

    fig.subplots_adjust(left=0.09, right=0.99, top=0.985, bottom=0.065, hspace=0.25, wspace=0.24)
    fig.savefig(SCRIPT_DIR / "supp_annotation_binning.pdf", format="pdf", bbox_inches="tight")
    fig.savefig(SCRIPT_DIR / "supp_annotation_binning.png", format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    plot()
    print(f"Wrote {SCRIPT_DIR / 'supp_annotation_binning.pdf'}")
    print(f"Wrote {SCRIPT_DIR / 'supp_annotation_binning.png'}")


if __name__ == "__main__":
    main()
