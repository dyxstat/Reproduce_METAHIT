#!/usr/bin/env python3
from pathlib import Path
import os

SCRIPT_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
import pandas as pd

PROJECT_DIR = SCRIPT_DIR.parents[1]
FONT_FAMILY = "Arial"
FONT_DIR = PROJECT_DIR / "conda_envs" / "metahit_env" / "fonts"
ARIAL_FILES = [FONT_DIR / name for name in ("arial.ttf", "arialbd.ttf", "arialbi.ttf", "ariali.ttf")]
for font_file in ARIAL_FILES:
    if not font_file.is_file():
        raise RuntimeError(f"Arial font file is missing: {font_file}")
    font_manager.fontManager.addfont(str(font_file))
font_manager.findfont(FONT_FAMILY, fallback_to_default=False)

DATASETS = [
    ("Human", "hg", PROJECT_DIR / "hg/10_MGE/results/genomad_output"),
    ("Pig", "pig", PROJECT_DIR / "pig/10_MGE/results/genomad_output"),
    ("Bovine", "bovine", PROJECT_DIR / "bovine/10_MGE/results/genomad_output"),
    ("Wastewater", "ww", PROJECT_DIR / "ww/10_MGE/results/genomad_output"),
    ("Mats", "mat", PROJECT_DIR / "mat/10_MGE/results/genomad_output"),
]

def find_summary(genomad_dir, kind):
    matches = sorted(genomad_dir.glob(f"**/*_{kind}_summary.tsv"))
    if not matches:
        raise FileNotFoundError(f"Missing geNomad {kind} summary under {genomad_dir}")
    if len(matches) > 1:
        summary_matches = [p for p in matches if p.parent.name.endswith("_summary")]
        if len(summary_matches) == 1:
            return summary_matches[0]
    return matches[0]


def count_summary(summary_tsv):
    df = pd.read_csv(summary_tsv, sep="\t")
    if "seq_name" not in df.columns:
        raise ValueError(f"Missing seq_name column in {summary_tsv}")
    names = df["seq_name"].astype(str)
    residual = int(names.str.startswith("residual|").sum())
    bin_specific = int(len(names) - residual)
    return residual, bin_specific


def load_counts():
    counts = {}
    missing = []
    for title, key, genomad_dir in DATASETS:
        try:
            virus_summary = find_summary(genomad_dir, "virus")
            plasmid_summary = find_summary(genomad_dir, "plasmid")
            counts[key] = {
                "title": title,
                "virus": count_summary(virus_summary),
                "plasmid": count_summary(plasmid_summary),
                "virus_summary": str(virus_summary),
                "plasmid_summary": str(plasmid_summary),
            }
        except (FileNotFoundError, ValueError) as exc:
            missing.append(str(exc))

    if missing:
        msg = "\n".join(missing)
        raise SystemExit(
            "Cannot generate all_pies.pdf/png because required geNomad summaries are missing:\n"
            f"{msg}\n"
            "Run Module 10/geNomad for the missing datasets, then rerun this script."
        )
    return counts


def make_plot(counts):
    matplotlib.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    virus_colors = ["#cce6ff", "#f5f5f5"]
    plasmid_colors = ["#ffcccc", "#f5f5f5"]

    fig, axes = plt.subplots(2, 5, figsize=(30, 12))

    for col, (title, key, _) in enumerate(DATASETS):
        ax_v = axes[0, col]
        ax_v.pie(
            counts[key]["virus"],
            colors=virus_colors,
            startangle=90,
            counterclock=True,
        )
        ax_v.axis("equal")
        ax_v.set_title(title, fontsize=36, fontname=FONT_FAMILY, fontweight="normal")

        ax_p = axes[1, col]
        ax_p.pie(
            counts[key]["plasmid"],
            colors=plasmid_colors,
            startangle=90,
            counterclock=True,
        )
        ax_p.axis("equal")

    legend_elements = [
        Patch(facecolor="#cce6ff", edgecolor="none", label="Viral contigs from residual assembly"),
        Patch(facecolor="#ffcccc", edgecolor="none", label="Plasmid contigs from residual assembly"),
        Patch(facecolor="#f5f5f5", edgecolor="none", label="contigs from bin-specific reassemblies"),
    ]

    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=3,
        fontsize=30,
        frameon=False,
    )

    plt.subplots_adjust(hspace=0.05, bottom=0.15)
    for suffix in ("pdf", "png"):
        output = SCRIPT_DIR / f"all_pies.{suffix}"
        fig.savefig(output, dpi=300, bbox_inches="tight")
        print(f"Saved {output}")
    plt.close(fig)


def main():
    counts = load_counts()
    for _, key, _ in DATASETS:
        item = counts[key]
        print(
            f"{item['title']}\t"
            f"virus_residual={item['virus'][0]}\tvirus_bin={item['virus'][1]}\t"
            f"plasmid_residual={item['plasmid'][0]}\tplasmid_bin={item['plasmid'][1]}"
        )
    make_plot(counts)


if __name__ == "__main__":
    main()
