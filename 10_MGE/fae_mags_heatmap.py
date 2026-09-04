#!/usr/bin/env python3
from pathlib import Path
import os

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import font_manager
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.sparse import load_npz

GTDB_SUMMARY = PROJECT_DIR / "hg/9_annotation/results/classify/gtdbtk.bac120.summary.tsv"
BIN_DIR = PROJECT_DIR / "hg/7_reassembly/results/reassembly_sg_hic/reassembled_bins"
CONTACT_MATRIX = PROJECT_DIR / "hg/10_MGE/results/10_MGE/mge_contact/mge_contact/denoised_contact_matrix_normcc.npz"
CONTIG_INFO = PROJECT_DIR / "hg/10_MGE/results/10_MGE/mge_contact/mge_contact/contig_info.csv"

OUT_PDF = SCRIPT_DIR / "fae_mags_heatmap.pdf"
OUT_PNG = SCRIPT_DIR / "fae_mags_heatmap.png"
OUT_TSV = SCRIPT_DIR / "fae_mags_heatmap_bins.tsv"
LABEL_FONT_SIZE = 42
RANDOM_SEED = 0
FONT_FAMILY = "Arial"
FONT_DIR = PROJECT_DIR / "conda_envs" / "metahit_env" / "fonts"
ARIAL_FILES = [FONT_DIR / name for name in ("arial.ttf", "arialbd.ttf", "arialbi.ttf", "ariali.ttf")]
for font_file in ARIAL_FILES:
    if not font_file.is_file():
        raise RuntimeError(f"Arial font file is missing: {font_file}")
    font_manager.fontManager.addfont(str(font_file))
font_manager.findfont(FONT_FAMILY, fallback_to_default=False)


def require_file(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input: {path}")


def fasta_ids(path):
    ids = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def extract_rank(lineage, prefix):
    for item in str(lineage).split(";"):
        if item.startswith(prefix):
            return item[len(prefix):]
    return ""


def label_from_bin(bin_name):
    core = bin_name.split(".", 2)
    if len(core) >= 2:
        return f"Bin{core[1]}"
    return bin_name


def load_faecalibacterium_bins():
    df = pd.read_csv(GTDB_SUMMARY, sep="\t")
    genus = df["classification"].map(lambda x: extract_rank(x, "g__"))
    fae = df[genus.eq("Faecalibacterium")].copy()
    if fae.empty:
        raise RuntimeError("No Faecalibacterium MAGs found in GTDB-Tk summary")

    rows = []
    for _, row in fae.iterrows():
        bin_name = str(row["user_genome"]).replace(".fa", "")
        fasta = BIN_DIR / f"{bin_name}.fa"
        require_file(fasta)
        local_contigs = fasta_ids(fasta)
        rows.append(
            {
                "bin": bin_name,
                "label": label_from_bin(bin_name),
                "species": extract_rank(row["classification"], "s__"),
                "contig_count": len(local_contigs),
                "fasta": fasta,
                "prefixed_contigs": [f"{bin_name}|{contig}" for contig in local_contigs],
            }
        )

    return rows


def make_heatmap(rows):
    contig_order = pd.read_csv(CONTIG_INFO)["name"].astype(str).tolist()
    contig_to_index = {contig: idx for idx, contig in enumerate(contig_order)}
    matrix = load_npz(CONTACT_MATRIX).tocsr()
    if matrix.shape[0] != len(contig_order) or matrix.shape[1] != len(contig_order):
        raise ValueError("Contact matrix dimensions do not match contig_info.csv")

    selected_rows = []
    for item in rows:
        present = [contig for contig in item["prefixed_contigs"] if contig in contig_to_index]
        missing = len(item["prefixed_contigs"]) - len(present)
        if not present:
            raise RuntimeError(f"{item['bin']} has no contigs present in the contact matrix")
        selected = dict(item)
        selected["present_contigs"] = present
        selected["missing_from_contact_matrix"] = missing
        selected_rows.append(selected)

    selected_rows.sort(key=lambda item: len(item["present_contigs"]), reverse=True)
    bin_matrices = []
    tick_locs = [0]
    labels = []
    output_rows = []

    rng = np.random.default_rng(RANDOM_SEED)
    for item in selected_rows:
        indices = [contig_to_index[contig] for contig in item["present_contigs"]]
        mat_bin = matrix[indices, :][:, indices].toarray()
        perm = rng.permutation(mat_bin.shape[0])
        mat_bin = mat_bin[perm, :][:, perm]
        bin_matrices.append(mat_bin)
        labels.append(item["label"])
        tick_locs.append(tick_locs[-1] + len(indices))
        output_rows.append(
            {
                "label": item["label"],
                "bin": item["bin"],
                "species": item["species"] if item["species"] else "unresolved",
                "total_contig_count": item["contig_count"],
                "contact_matrix_contig_count": len(indices),
                "missing_from_contact_matrix": item["missing_from_contact_matrix"],
                "fasta": str(item["fasta"]),
            }
        )

    pd.DataFrame(output_rows).to_csv(OUT_TSV, sep="\t", index=False)

    submatrix = np.block(
        [
            [
                bin_matrices[i]
                if i == j
                else np.zeros((bin_matrices[i].shape[0], bin_matrices[j].shape[0]))
                for j in range(len(bin_matrices))
            ]
            for i in range(len(bin_matrices))
        ]
    )
    np.fill_diagonal(submatrix, 0)
    submatrix = np.log(submatrix + 0.01)

    matplotlib.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, ax = plt.subplots(figsize=(18, 16))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.25)
    heat = sns.heatmap(
        submatrix,
        square=True,
        cmap="rocket",
        ax=ax,
        cbar=True,
        cbar_ax=cax,
        linewidths=0,
        vmin=-4,
        vmax=8,
    )
    cbar = heat.collections[0].colorbar
    cbar.ax.tick_params(labelsize=LABEL_FONT_SIZE)

    for loc in tick_locs:
        ax.hlines(loc, *ax.get_xlim(), color="grey", linewidth=0.5, linestyle="-.")
        ax.vlines(loc, *ax.get_ylim(), color="grey", linewidth=0.5, linestyle="-.")

    midpoints = [(tick_locs[i] + tick_locs[i + 1]) / 2 for i in range(len(labels))]
    ax.set_yticks(midpoints)
    ax.set_yticklabels(labels, rotation=0, fontsize=LABEL_FONT_SIZE)
    ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax.set_xlabel("")
    ax.set_ylabel("")

    plt.tight_layout(pad=1.5)
    fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight", pad_inches=0.25)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)


def main():
    for path in [GTDB_SUMMARY, CONTACT_MATRIX, CONTIG_INFO]:
        require_file(path)
    rows = load_faecalibacterium_bins()
    make_heatmap(rows)
    print(f"Saved {OUT_TSV}")
    print(f"Saved {OUT_PDF}")
    print(f"Saved {OUT_PNG}")
    print(pd.read_csv(OUT_TSV, sep="\t").to_string(index=False))


if __name__ == "__main__":
    main()
