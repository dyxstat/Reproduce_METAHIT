#!/usr/bin/env python3
import os
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
MPL_DIR = SCRIPT_DIR / ".matplotlib"
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager


QUALITY_MIN_COMPLETENESS = 50.0
QUALITY_MAX_CONTAMINATION = 10.0

CSV_FILE = SCRIPT_DIR / "phylum_stacked.csv"
PDF_FILE = SCRIPT_DIR / "phylum_stacked.pdf"
PNG_FILE = SCRIPT_DIR / "phylum_stacked.png"

CURRENT_DATASETS = [
    {
        "dataset": "human gut",
        "quality": ROOT / "hg/7_reassembly/results/reassembly_sg_hic/reassembled_bins.checkm2/quality_report.tsv",
        "tax_dir": ROOT / "hg/9_annotation/results/classify",
    },
    {
        "dataset": "sheep_gut",
        "quality": ROOT / "sheep/6_binning/results/metahit/metahit_50_10_bins.stats",
        "tax_dir": ROOT / "sheep/9_annotation/results/classify",
    },
    {
        "dataset": "pig_gut",
        "quality": ROOT / "pig/7_reassembly/results/reassembly_sg_hic/reassembled_bins.checkm2/quality_report.tsv",
        "tax_dir": ROOT / "pig/9_annotation/results/classify",
    },
    {
        "dataset": "cow_rumen",
        "quality": ROOT / "cow/6_binning/results/metahit/work_files/binsO.checkm2/quality_report.tsv",
        "tax_dir": ROOT / "cow/9_annotation/results/classify",
    },
    {
        "dataset": "bovine_skin",
        "quality": ROOT / "bovine/7_reassembly/results/reassembly_sg_hic/reassembled_bins.checkm2/quality_report.tsv",
        "tax_dir": ROOT / "bovine/9_annotation/results/classify",
    },
    {
        "dataset": "wastewater",
        "quality": ROOT / "ww/7_reassembly/results/reassembly_sg_hic/reassembled_bins.checkm2/quality_report.tsv",
        "tax_dir": ROOT / "ww/9_annotation/results/classify",
    },
    {
        "dataset": "hydrothermal_mats",
        "quality": ROOT / "mat/7_reassembly/results/reassembly_sg_hic/reassembled_bins.checkm2/quality_report.tsv",
        "tax_dir": ROOT / "mat/9_annotation/results/classify",
    },
]

DATASET_ORDER = [
    "human gut",
    "sheep_gut",
    "pig_gut",
    "cow_rumen",
    "bovine_skin",
    "wastewater",
    "hydrothermal_mats",
]

RENAME_MAP = {
    "human gut": "Human",
    "sheep_gut": "Sheep",
    "pig_gut": "Pig",
    "cow_rumen": "Cow",
    "bovine_skin": "Bovine",
    "wastewater": "Wastewater",
    "hydrothermal_mats": "Mats",
}

FONT_FAMILY = "Arial"
FONT_DIR = ROOT / "conda_envs" / "metahit_env" / "fonts"
ARIAL_FILES = [FONT_DIR / name for name in ("arial.ttf", "arialbd.ttf", "arialbi.ttf", "ariali.ttf")]
for font_file in ARIAL_FILES:
    if not font_file.is_file():
        raise RuntimeError(f"Arial font file is missing: {font_file}")
    font_manager.fontManager.addfont(str(font_file))
font_manager.findfont(FONT_FAMILY, fallback_to_default=False)


def extract_phylum(classification):
    if pd.isna(classification):
        return "Unclassified"
    for field in str(classification).split(";"):
        if field.startswith("p__"):
            return field[3:] or "Unclassified"
    return "Unclassified"


def read_taxonomy(tax_dir):
    frames = []
    for domain in ("bac120", "ar53"):
        path = tax_dir / f"gtdbtk.{domain}.summary.tsv"
        if path.exists():
            frames.append(pd.read_csv(path, sep="\t")[["user_genome", "classification"]])
    if not frames:
        raise FileNotFoundError(f"No GTDB-Tk summary files found in {tax_dir}")
    tax = pd.concat(frames, ignore_index=True)
    tax["user_genome"] = tax["user_genome"].astype(str)
    return tax.drop_duplicates(subset=["user_genome"], keep="first")


def read_quality(quality_path):
    if not quality_path.exists():
        raise FileNotFoundError(f"Missing CheckM2 quality report: {quality_path}")
    quality = pd.read_csv(quality_path, sep="\t")
    column_lookup = {column.lower(): column for column in quality.columns}
    needed = {"name", "completeness", "contamination"}
    if "name" not in column_lookup and "bin" in column_lookup:
        column_lookup["name"] = column_lookup["bin"]
    missing = needed - set(column_lookup)
    if missing:
        raise ValueError(f"{quality_path} is missing columns: {sorted(missing)}")
    quality = quality.rename(
        columns={
            column_lookup["name"]: "Name",
            column_lookup["completeness"]: "Completeness",
            column_lookup["contamination"]: "Contamination",
        }
    )
    quality = quality[
        (quality["Completeness"] >= QUALITY_MIN_COMPLETENESS)
        & (quality["Contamination"] < QUALITY_MAX_CONTAMINATION)
    ].copy()
    quality["Name"] = quality["Name"].astype(str)
    return quality


def percentages_from_outputs(dataset):
    quality = read_quality(dataset["quality"])
    tax = read_taxonomy(dataset["tax_dir"])
    merged = quality.merge(tax, left_on="Name", right_on="user_genome", how="left")

    missing_tax = merged.loc[merged["classification"].isna(), "Name"].tolist()
    if missing_tax:
        preview = ", ".join(missing_tax[:10])
        raise RuntimeError(
            f"Missing GTDB-Tk taxonomy for {dataset['dataset']} bins: {preview}"
            + (" ..." if len(missing_tax) > 10 else "")
        )

    if merged.empty:
        return {}
    phylum = merged["classification"].map(extract_phylum)
    return (phylum.value_counts() / len(phylum) * 100.0).to_dict()


def build_updated_table():
    current_values = {
        dataset["dataset"]: percentages_from_outputs(dataset) for dataset in CURRENT_DATASETS
    }
    primary_phyla = [
        "Bacteroidota", "Bacillota_A", "Pseudomonadota", "Campylobacterota",
        "Spirochaetota", "Bacillota", "Actinomycetota", "Desulfobacterota",
    ]
    other_phyla = sorted(
        set().union(*(set(values) for values in current_values.values())) - set(primary_phyla)
    )
    columns = ["dataset", *primary_phyla, *other_phyla]
    current_rows = []
    for dataset in CURRENT_DATASETS:
        row = {col: 0.0 for col in columns}
        row["dataset"] = dataset["dataset"]
        row.update(current_values[dataset["dataset"]])
        current_rows.append(row)

    combined = pd.DataFrame(current_rows, columns=columns)
    combined = combined.set_index("dataset").loc[DATASET_ORDER].reset_index()
    return combined


def make_plot(df):
    plt.rcParams["font.family"] = FONT_FAMILY
    plt.rcParams["font.size"] = 50
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    df = df.copy()
    df["dataset"] = df["dataset"].replace(RENAME_MAP)

    left_cols = df.columns[1:9]

    df_plot = df.copy()
    df_plot["others"] = 100 - df[left_cols].sum(axis=1)
    df_plot["others"] = df_plot["others"].clip(lower=0)
    plot_df = df_plot[["dataset"] + list(left_cols) + ["others"]]

    custom_order = ["Human", "Sheep", "Pig", "Cow", "Bovine", "Wastewater", "Mats"]
    plot_df = plot_df.set_index("dataset").loc[custom_order]

    palette_named = [
        "#0072B2",
        "#D55E00",
        "#009E73",
        "#CC79A7",
        "#E69F00",
        "#56B4E9",
        "#F0E442",
        "#6F4C9B",
        "#1B9E77",
        "#E31A1C",
        "#A6CEE3",
        "#33A02C",
    ]
    colors_lookup = {}
    for i, col in enumerate(plot_df.columns):
        if col == "others":
            colors_lookup[col] = "#9A9A9A"
        else:
            colors_lookup[col] = palette_named[i % len(palette_named)]
    color_sequence = [colors_lookup[c] for c in plot_df.columns]

    ax = plot_df.plot(
        kind="bar",
        stacked=True,
        figsize=(36, 20),
        width=0.6,
        color=color_sequence,
    )

    ax.set_xlabel("Environment", fontsize=50)
    ax.set_ylabel("Percentage", fontsize=50)
    ax.tick_params(axis="x", labelsize=50)
    ax.tick_params(axis="y", labelsize=50)
    ax.set_ylim(0, 100)

    plt.legend(
        bbox_to_anchor=(0.5, -0.08),
        loc="upper center",
        ncol=(len(plot_df.columns) + 1) // 2,
        frameon=False,
        fontsize=40,
    )

    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(PDF_FILE, bbox_inches="tight")
    plt.savefig(PNG_FILE, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    df = build_updated_table()
    df.to_csv(CSV_FILE, index=False, float_format="%.10g")
    make_plot(df)

    print(f"Wrote updated CSV: {CSV_FILE}")
    print(f"Wrote PDF: {PDF_FILE}")
    print(f"Wrote PNG: {PNG_FILE}")


if __name__ == "__main__":
    main()
