#!/usr/bin/env python3
import os
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".matplotlib"))

import matplotlib.pyplot as plt
from matplotlib import font_manager
import seaborn as sns

PROJECT_DIR = SCRIPT_DIR.parents[1]
FONT_FAMILY = "Arial"
FONT_DIR = PROJECT_DIR / "conda_envs" / "metahit_env" / "fonts"
ARIAL_FILES = [FONT_DIR / name for name in ("arial.ttf", "arialbd.ttf", "arialbi.ttf", "ariali.ttf")]
for font_file in ARIAL_FILES:
    if not font_file.is_file():
        raise RuntimeError(f"Arial font file is missing: {font_file}")
    font_manager.fontManager.addfont(str(font_file))
font_manager.findfont(FONT_FAMILY, fallback_to_default=False)


DATASET_RENAME = {
    "human gut": "Human",
    "sheep_gut": "Sheep",
    "pig_gut": "Pig",
    "cow_rumen": "Cow",
    "bovine_skin": "Bovine",
    "wastewater": "Wastewater",
    "hydrothermal_mats": "Mats",
}

SAMPLES = ["Human", "Sheep", "Pig", "Cow", "Bovine", "Wastewater", "Mats"]


def bray_curtis(x, y):
    denominator = np.sum(x + y)
    if denominator == 0:
        return 0.0
    return np.sum(np.abs(x - y)) / denominator


def load_phylum_table(csv_file):
    df = pd.read_csv(csv_file)
    if "dataset" not in df.columns:
        raise ValueError("Input CSV must contain a 'dataset' column.")
    df["dataset"] = df["dataset"].replace(DATASET_RENAME)
    df = df.set_index("dataset")
    missing = [sample for sample in SAMPLES if sample not in df.index]
    if missing:
        raise ValueError(f"Missing dataset(s) in {csv_file}: {', '.join(missing)}")
    return df.loc[SAMPLES].astype(float)


def compute_matrix(phylum_df):
    matrix = np.zeros((len(SAMPLES), len(SAMPLES)), dtype=float)
    values = phylum_df.to_numpy()
    for i in range(len(SAMPLES)):
        for j in range(len(SAMPLES)):
            matrix[i, j] = bray_curtis(values[i], values[j])
    return pd.DataFrame(matrix, index=SAMPLES, columns=SAMPLES)


def plot_dissimilarity(df, pdf_file, png_file):
    plt.rcParams["font.family"] = FONT_FAMILY

    mask = np.triu(np.ones_like(df, dtype=bool), k=0)
    plt.figure(figsize=(10, 8))

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(color="white")

    ax = sns.heatmap(
        df,
        annot=False,
        cmap=cmap,
        mask=mask,
        square=True,
        linewidths=0.5,
        vmin=0,
        vmax=1,
        cbar_kws={"label": "", "shrink": 0.8},
    )

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=25)

    plt.xlabel("")
    plt.ylabel("")
    plt.xticks(rotation=45, ha="right", fontsize=25)
    plt.yticks(rotation=0, fontsize=25)
    plt.tight_layout()
    plt.savefig(pdf_file, dpi=300, bbox_inches="tight")
    plt.savefig(png_file, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    csv_file = SCRIPT_DIR / "phylum_stacked.csv"
    matrix_file = SCRIPT_DIR / "dissimilarity_matrix.csv"
    pdf_file = SCRIPT_DIR / "dissimilarity.pdf"
    png_file = SCRIPT_DIR / "dissimilarity.png"

    phylum_df = load_phylum_table(csv_file)
    dissimilarity_df = compute_matrix(phylum_df)
    dissimilarity_df.round(3).to_csv(matrix_file)
    plot_dissimilarity(dissimilarity_df, pdf_file, png_file)

    print("Dissimilarity matrix:")
    print(dissimilarity_df.round(3).to_string())
    print(f"Wrote matrix: {matrix_file}")
    print(f"Wrote PDF: {pdf_file}")
    print(f"Wrote PNG: {png_file}")


if __name__ == "__main__":
    main()
