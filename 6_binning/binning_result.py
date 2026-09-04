#!/usr/bin/env python3
"""Plot Hi-C binner comparison counts from the seven stored result sets.

All counts are calculated directly from the per-dataset ``results`` statistics
files.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[1]
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

METHODS = ["METAHICT", "bin3C", "MetaCC", "ImputeCC"]
DATASET_ORDER = [
    ("hg", "Human Gut"),
    ("sheep", "Sheep Gut"),
    ("pig", "Pig Gut"),
    ("cow", "Cow Rumen"),
    ("bovine", "Bovine Skin"),
    ("ww", "Wastewater"),
    ("mat", "Hydrothermal Mats"),
]
METHOD_STAT_FILES = [
    ("METAHICT", "metahit_50_10_bins.stats"),
    ("bin3C", "work_files/binsC.stats"),
    ("MetaCC", "work_files/binsA.stats"),
    ("ImputeCC", "work_files/binsB.stats"),
]
RESULTS_DIRECTORIES = {
    "hg": PROJECT_DIR / "hg" / "6_binning" / "results" / "metahit",
    "sheep": PROJECT_DIR / "sheep" / "6_binning" / "results" / "metahit",
    "pig": PROJECT_DIR / "pig" / "6_binning" / "results" / "metahit",
    "cow": PROJECT_DIR / "cow" / "6_binning" / "results" / "metahit",
    "bovine": PROJECT_DIR / "bovine" / "6_binning" / "results" / "metahit",
    "ww": PROJECT_DIR / "ww" / "6_binning" / "results" / "metahit",
    "mat": PROJECT_DIR / "mat" / "6_binning" / "results" / "metahit",
}
COMP_THRESHOLDS = [50, 70, 90]
CONT_THRESHOLDS = [5, 10]
BLUE_COLORS = {50: "#C7D2E3", 70: "#8198C0", 90: "#4F74B7"}
RED_COLORS = {50: "#F0C2C1", 70: "#E3918F", 90: "#DD514A"}


def tick_interval(max_value: int) -> int:
    """Return x-axis tick spacing matching the old reference style."""
    if max_value >= 500:
        return 200
    if max_value >= 150:
        return 50
    return 20


def method_counts(values: list[int]) -> np.ndarray:
    return np.array(values)


def read_stats(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing stats file: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        completeness_field = next((name for name in fieldnames if name.lower() == "completeness"), None)
        contamination_field = next((name for name in fieldnames if name.lower() == "contamination"), None)
        missing = [
            label
            for label, field in (
                ("completeness", completeness_field),
                ("contamination", contamination_field),
            )
            if field is None
        ]
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
        return [
            {
                "completeness": row[completeness_field],
                "contamination": row[contamination_field],
            }
            for row in reader
        ]


def count_bins(rows: list[dict[str, str]], min_completeness: int, max_contamination: int) -> int:
    return sum(
        float(row["completeness"]) >= min_completeness
        and float(row["contamination"]) < max_contamination
        for row in rows
    )


def dataset_counts(dataset_key: str) -> dict[int, dict[int, np.ndarray]]:
    """Return threshold counts for one dataset.

    The METAHICT integration script writes its three input binners as:
    binsA = MetaCC, binsB = ImputeCC, and binsC = bin3C.  The order below is
    kept consistent with METHODS for plotting: METAHICT, bin3C, MetaCC,
    ImputeCC.
    """
    stats_dir = RESULTS_DIRECTORIES[dataset_key]
    rows_by_method = [
        read_stats(stats_dir / relative_path)
        for method_label, relative_path in METHOD_STAT_FILES
    ]
    counts: dict[int, dict[int, np.ndarray]] = {}
    for cont in CONT_THRESHOLDS:
        counts[cont] = {}
        for comp in COMP_THRESHOLDS:
            values = [count_bins(rows, comp, cont) for rows in rows_by_method]
            counts[cont][comp] = method_counts(values)
    return counts


DATASETS = [
    (dataset_label, dataset_counts(dataset_key))
    for dataset_key, dataset_label in DATASET_ORDER
]


def write_label_map(path: Path) -> None:
    with path.open("w") as handle:
        for idx, (dataset_label, _data) in enumerate(DATASETS):
            handle.write(f"{chr(ord('A') + idx)}: {dataset_label}\n")


def plot() -> None:
    plt.rcParams["font.family"] = FONT_FAMILY
    fig, axes = plt.subplots(nrows=len(DATASETS), ncols=2, figsize=(14, 22), sharey=True)

    y = np.arange(len(METHODS))
    method_labels = METHODS

    for row_idx, (_dataset_label, data) in enumerate(DATASETS):
        for col_idx, cont in enumerate(CONT_THRESHOLDS):
            ax = axes[row_idx, col_idx]
            colors = BLUE_COLORS if cont == 5 else RED_COLORS
            # Match the reference Plotly figure more closely: each subplot
            # should autoscale from its own values.  The previous matplotlib
            # version used the larger of the two columns in the row plus 12%
            # padding, which stretched the left HG panel to a 140 tick.
            col_max = int(data[cont][50].max())
            x_max = max(1, col_max * 1.03)
            tick_step = tick_interval(col_max)

            for comp in [50, 70, 90]:
                ax.barh(
                    y,
                    data[cont][comp],
                    height=0.68,
                    color=colors[comp],
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
                    chr(ord("A") + row_idx),
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

    fig.savefig(SCRIPT_DIR / "binning_result.pdf", format="pdf", bbox_inches="tight")
    fig.savefig(SCRIPT_DIR / "binning_results.png", format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    write_label_map(SCRIPT_DIR / "binning_results.txt")
    plot()
    print(f"Wrote {SCRIPT_DIR / 'binning_result.pdf'}")
    print(f"Wrote {SCRIPT_DIR / 'binning_results.png'}")
    print(f"Wrote {SCRIPT_DIR / 'binning_results.txt'}")


if __name__ == "__main__":
    main()
