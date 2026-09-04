#!/usr/bin/env python3
"""Plot direct CheckM2 threshold counts for Hi-C and non-Hi-C binners."""

from __future__ import annotations

import csv
import os
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[1]
BASELINE_DIR = PROJECT_DIR / "baseline_results"
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".matplotlib"))

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MultipleLocator
import numpy as np


FONT_FAMILY = "Arial"
FONT_DIR = PROJECT_DIR / "conda_envs" / "metahit_env" / "fonts"
ARIAL_FILES = [FONT_DIR / name for name in ("arial.ttf", "arialbd.ttf", "arialbi.ttf", "ariali.ttf")]
for font_file in ARIAL_FILES:
    if not font_file.is_file():
        raise RuntimeError(f"Arial font file is missing: {font_file}")
    font_manager.fontManager.addfont(str(font_file))
font_manager.findfont(FONT_FAMILY, fallback_to_default=False)

DATASETS = [
    ("hg", "Human Gut"),
    ("sheep", "Sheep Gut"),
    ("pig", "Pig Gut"),
    ("cow", "Cow Rumen"),
    ("bovine", "Bovine Skin"),
    ("ww", "Wastewater"),
    ("mat", "Hydrothermal Mats"),
]

METHODS = [
    ("bin3c", "bin3C", ("binning_refiner", "binning_refiner", "output", "checkm2", "bin3c", "quality_report.tsv")),
    ("metacc", "MetaCC", ("binning_refiner", "binning_refiner", "output", "checkm2", "metacc", "quality_report.tsv")),
    ("imputecc", "ImputeCC", ("binning_refiner", "binning_refiner", "output", "checkm2", "imputecc", "quality_report.tsv")),
    ("metabat2", "MetaBAT2", ("metabat2", "output", "checkm2", "quality_report.tsv")),
    ("concoct", "CONCOCT", ("concoct", "output", "checkm2", "quality_report.tsv")),
]

COMP_THRESHOLDS = [50, 70, 90]
CONT_THRESHOLDS = [5, 10]

COLORS = {
    5: {50: "#C7D2E3", 70: "#8198C0", 90: "#4F74B7"},
    10: {50: "#F0C2C1", 70: "#E3918F", 90: "#DD514A"},
}


def tick_interval(max_value: int) -> int:
    """Return x-axis tick spacing matching update_plot/binning_result.py."""
    if max_value >= 500:
        return 200
    if max_value >= 150:
        return 50
    return 20

def read_checkm2_report(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing CheckM2 report: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Completeness", "Contamination"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
        return list(reader)


def count_bins(rows: list[dict[str, str]], min_completeness: int, max_contamination: int) -> int:
    return sum(
        float(row["Completeness"]) >= min_completeness
        and float(row["Contamination"]) < max_contamination
        for row in rows
    )


def collect_counts() -> dict[tuple[str, str, int, int], int]:
    counts: dict[tuple[str, str, int, int], int] = {}
    for dataset_key, _dataset_label in DATASETS:
        for method_key, _method_label, report_parts in METHODS:
            report = BASELINE_DIR / dataset_key / Path(*report_parts)
            rows = read_checkm2_report(report)
            for cont in CONT_THRESHOLDS:
                for comp in COMP_THRESHOLDS:
                    counts[(dataset_key, method_key, cont, comp)] = count_bins(rows, comp, cont)
    return counts


def write_label_map(output_txt: Path) -> None:
    with output_txt.open("w") as handle:
        for idx, (_dataset_key, dataset_label) in enumerate(DATASETS):
            letter = chr(ord("A") + idx)
            handle.write(f"{letter}: {dataset_label}\n")


def plot_counts(
    counts: dict[tuple[str, str, int, int], int],
    output_pdf: Path,
    output_png: Path,
) -> None:
    plt.rcParams.update({
        "font.family": [FONT_FAMILY],
        "font.sans-serif": [FONT_FAMILY],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, axes = plt.subplots(
        nrows=len(DATASETS),
        ncols=2,
        figsize=(14, 22),
        sharey=True,
    )

    y = np.arange(len(METHODS))
    method_labels = [label for _key, label, _parts in METHODS]
    letters = [chr(ord("A") + i) for i in range(len(DATASETS))]

    for row_idx, ((dataset_key, dataset_label), letter) in enumerate(zip(DATASETS, letters)):
        for col_idx, cont in enumerate(CONT_THRESHOLDS):
            ax = axes[row_idx, col_idx]
            col_max = max(
                counts[(dataset_key, method_key, cont, 50)]
                for method_key, _method_label, _parts in METHODS
            )
            x_max = max(1, col_max * 1.03)
            tick_step = tick_interval(col_max)

            for comp in COMP_THRESHOLDS:
                values = [
                    counts[(dataset_key, method_key, cont, comp)]
                    for method_key, _method_label, _parts in METHODS
                ]
                ax.barh(
                    y,
                    values,
                    height=0.68,
                    color=COLORS[cont][comp],
                    edgecolor="white",
                    linewidth=0.7,
                    label=f"Comp ≥ {comp}%",
                    zorder=comp,
                )

            ax.set_xlim(0, x_max)
            ax.xaxis.set_major_locator(MultipleLocator(tick_step))
            ax.set_yticks(y)
            ax.set_yticklabels(method_labels, fontsize=13)
            ax.set_ylim(len(METHODS) - 0.5, -0.5)
            ax.set_xlabel("Number of bins", fontsize=13)
            ax.tick_params(axis="x", labelsize=12)
            ax.grid(axis="x", color="#E2E2E2", linewidth=0.8)
            ax.set_axisbelow(True)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            if col_idx == 0:
                ax.text(
                    -0.28,
                    -0.16,
                    letter,
                    transform=ax.get_yaxis_transform(),
                    fontsize=34,
                    fontfamily=FONT_FAMILY,
                    fontweight="normal",
                    va="top",
                    ha="right",
                )

    fig.subplots_adjust(left=0.24, right=0.98, top=0.965, bottom=0.065, hspace=0.56, wspace=0.22)
    left_bbox = axes[0, 0].get_position()
    right_bbox = axes[0, 1].get_position()
    left_legend_x = (left_bbox.x0 + left_bbox.x1) / 2
    right_legend_x = (right_bbox.x0 + right_bbox.x1) / 2

    handles_left, labels_left = axes[0, 0].get_legend_handles_labels()
    handles_right, labels_right = axes[0, 1].get_legend_handles_labels()
    fig.legend(
        handles_left,
        labels_left,
        title="Cont < 5%",
        loc="lower center",
        bbox_to_anchor=(left_legend_x, 0.011),
        ncol=3,
        frameon=False,
        fontsize=12,
        title_fontsize=14,
    )
    fig.legend(
        handles_right,
        labels_right,
        title="Cont < 10%",
        loc="lower center",
        bbox_to_anchor=(right_legend_x, 0.011),
        ncol=3,
        frameon=False,
        fontsize=12,
        title_fontsize=14,
    )

    fig.savefig(output_pdf, format="pdf", bbox_inches="tight")
    fig.savefig(output_png, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    counts = collect_counts()
    output_pdf = SCRIPT_DIR / "non_hic_binners.pdf"
    output_png = SCRIPT_DIR / "non_hic_binners.png"
    output_txt = SCRIPT_DIR / "non_hic_binners.txt"
    write_label_map(output_txt)
    plot_counts(counts, output_pdf, output_png)
    print(f"Wrote {output_txt}")
    print(f"Wrote {output_pdf}")
    print(f"Wrote {output_png}")


if __name__ == "__main__":
    main()
