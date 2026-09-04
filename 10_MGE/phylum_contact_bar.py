#!/usr/bin/env python3
from pathlib import Path
import os

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd
from scipy.sparse import load_npz

HG_DIR = PROJECT_DIR / "hg"
GTDB_SUMMARY = HG_DIR / "9_annotation/results/classify/gtdbtk.bac120.summary.tsv"
CONTIG_INFO = HG_DIR / "10_MGE/results/10_MGE/mge_contact/mge_contact/contig_info.csv"
NORMALIZED_CONTACT = HG_DIR / "10_MGE/results/10_MGE/mge_contact/mge_contact/denoised_contact_matrix_normcc.npz"
RAW_CONTACT = HG_DIR / "10_MGE/results/10_MGE/mge_contact/mge_contact/Raw_contact_matrix.npz"
VIRAL_FASTA = HG_DIR / "10_MGE/results/10_MGE/mge/mge_reports/virus_no_provirus.fna"
PLASMID_FASTA = (
    HG_DIR
    / "10_MGE/results/10_MGE/mge/genomad_output/combined_contigs.unique_summary/combined_contigs.unique_plasmid.fna"
)
ASSOCIATION_TABLE = HG_DIR / "10_MGE/results/10_MGE/mge/candidate_mge_mag_associations_zscore_filtered.tsv"

MIN_RAW_CONTACTS = 2
ZSCORE_THRESHOLD = 0.5

OUT_TSV = SCRIPT_DIR / "phylum_contact_bar.tsv"
OUT_PDF = SCRIPT_DIR / "phylum_contact_bar.pdf"
OUT_PNG = SCRIPT_DIR / "phylum_contact_bar.png"
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
    ids = set()
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def load_bin_phyla(path):
    df = pd.read_csv(path, sep="\t")
    bin_to_phylum = {}
    for _, row in df.iterrows():
        bin_name = str(row["user_genome"]).replace(".fa", "")
        lineage = str(row["classification"]).split(";")
        phylum = next((item[3:] for item in lineage if item.startswith("p__")), None)
        if phylum:
            bin_to_phylum[bin_name] = phylum
            parts = bin_name.split(".")
            if len(parts) >= 2 and parts[0] == "bin" and parts[1].isdigit():
                bin_to_phylum[f"bin{parts[1]}"] = phylum
    return bin_to_phylum


def host_mag_from_contig(contig):
    if not contig.startswith("bin"):
        return None
    return contig.split("|", 1)[0]


def aggregate_scores(contigs, normalized_mat, raw_mat, viral_ids, plasmid_ids, bin_to_phylum):
    mge_type_by_id = {}
    for contig in viral_ids:
        mge_type_by_id[contig] = "Viral"
    for contig in plasmid_ids:
        mge_type_by_id.setdefault(contig, "Plasmid")

    associations = {}

    def add_normalized(mge_contig, host_contig, value):
        mge_type = mge_type_by_id.get(mge_contig)
        host_mag = host_mag_from_contig(host_contig)
        if not mge_type or not host_mag:
            return
        phylum = bin_to_phylum.get(host_mag)
        if not phylum:
            return
        key = (mge_type, mge_contig, host_mag)
        item = associations.setdefault(
            key,
            {
                "mge_type": mge_type,
                "mge_contig": mge_contig,
                "host_mag": host_mag,
                "phylum": phylum,
                "normalized_score": 0.0,
                "raw_support": 0.0,
            },
        )
        item["normalized_score"] += float(value)

    for r, c, value in zip(normalized_mat.row, normalized_mat.col, normalized_mat.data):
        if r >= c or float(value) <= 0:
            continue
        a, b = contigs[r], contigs[c]
        add_normalized(a, b, value)
        add_normalized(b, a, value)

    def add_raw(mge_contig, host_contig, value):
        mge_type = mge_type_by_id.get(mge_contig)
        host_mag = host_mag_from_contig(host_contig)
        if not mge_type or not host_mag:
            return
        key = (mge_type, mge_contig, host_mag)
        if key in associations:
            associations[key]["raw_support"] += float(value)

    for r, c, value in zip(raw_mat.row, raw_mat.col, raw_mat.data):
        if r >= c or float(value) <= 0:
            continue
        a, b = contigs[r], contigs[c]
        add_raw(a, b, value)
        add_raw(b, a, value)

    df = pd.DataFrame(associations.values())
    if df.empty:
        return df

    df["z_score"] = 0.0
    for mge_type, idx in df.groupby("mge_type").groups.items():
        scores = df.loc[idx, "normalized_score"]
        std = scores.std(ddof=0)
        if std and not pd.isna(std):
            df.loc[idx, "z_score"] = (scores - scores.mean()) / std

    return df[(df["raw_support"] >= MIN_RAW_CONTACTS) & (df["z_score"] > ZSCORE_THRESHOLD)].copy()


def summarize_by_phylum(associations):
    if associations.empty:
        return pd.DataFrame(columns=["Phylum", "Viral", "Plasmid", "Bins"])

    totals = associations.groupby(["phylum", "mge_type"])["mge_contig"].nunique().unstack(fill_value=0)
    for col in ["Viral", "Plasmid"]:
        if col not in totals.columns:
            totals[col] = 0
    totals["Total"] = totals["Viral"] + totals["Plasmid"]
    top_phyla = totals.sort_values("Total", ascending=False).head(4).index.tolist()

    rows = []
    for phylum in top_phyla:
        sub = associations[associations["phylum"] == phylum]
        rows.append(
            {
                "Phylum": phylum,
                "Viral": int(sub.loc[sub["mge_type"] == "Viral", "mge_contig"].nunique()),
                "Plasmid": int(sub.loc[sub["mge_type"] == "Plasmid", "mge_contig"].nunique()),
                "Bins": ",".join(sorted(sub["host_mag"].unique())),
            }
        )
    return pd.DataFrame(rows)


def load_associations_table(path, bin_to_phylum):
    df = pd.read_csv(path, sep="\t")
    required = {"mge_contig", "mge_type", "host_mag"}
    if not required.issubset(df.columns):
        return None

    associations = df.copy()
    associations["mge_type"] = (
        associations["mge_type"]
        .astype(str)
        .str.lower()
        .map({"virus": "Viral", "viral": "Viral", "plasmid": "Plasmid"})
    )
    associations = associations[associations["mge_type"].isin(["Viral", "Plasmid"])].copy()
    associations["phylum"] = associations["host_mag"].astype(str).map(bin_to_phylum)
    associations = associations[associations["phylum"].notna()].copy()
    return associations


def make_plot(df):
    matplotlib.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    plot_df = df.copy()
    plot_df["Total"] = plot_df["Viral"] + plot_df["Plasmid"]
    plot_df = plot_df.sort_values("Total", ascending=False)

    fig, ax = plt.subplots(figsize=(20, 12))
    x = range(len(plot_df))
    plasmid = plot_df["Plasmid"].astype(int).tolist()
    viral = plot_df["Viral"].astype(int).tolist()

    p1 = ax.bar(x, plasmid, color="#ffcccc", width=0.6, label="Plasmid contigs")
    p2 = ax.bar(x, viral, bottom=plasmid, color="#cce6ff", width=0.6, label="Viral contigs")

    ax.set_xticks(list(x))
    ax.set_xticklabels(plot_df["Phylum"].tolist(), fontsize=36, ha="center")
    ax.set_ylabel("Number of contigs", fontsize=40)
    ax.tick_params(axis="y", labelsize=36)

    fig.legend(
        handles=[p1, p2],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        ncol=2,
        fontsize=36,
        frameon=False,
    )
    plt.tight_layout(rect=[0, 0.12, 1, 1])
    fig.savefig(OUT_PDF, dpi=300)
    fig.savefig(OUT_PNG, dpi=300)
    plt.close(fig)


def main():
    for path in [GTDB_SUMMARY]:
        require_file(path)

    bin_to_phylum = load_bin_phyla(GTDB_SUMMARY)
    associations = None
    if ASSOCIATION_TABLE.is_file() and ASSOCIATION_TABLE.stat().st_size > 0:
        associations = load_associations_table(ASSOCIATION_TABLE, bin_to_phylum)

    if associations is None:
        for path in [CONTIG_INFO, NORMALIZED_CONTACT, RAW_CONTACT, VIRAL_FASTA, PLASMID_FASTA]:
            require_file(path)
        contigs = pd.read_csv(CONTIG_INFO)["name"].astype(str).tolist()
        normalized_mat = load_npz(NORMALIZED_CONTACT).tocoo()
        raw_mat = load_npz(RAW_CONTACT).tocoo()
        if normalized_mat.shape[0] != len(contigs) or raw_mat.shape[0] != len(contigs):
            raise ValueError("Contact matrix dimensions do not match contig_info.csv")

        associations = aggregate_scores(
            contigs,
            normalized_mat,
            raw_mat,
            fasta_ids(VIRAL_FASTA),
            fasta_ids(PLASMID_FASTA),
            bin_to_phylum,
        )
    summary = summarize_by_phylum(associations)
    summary.to_csv(OUT_TSV, sep="\t", index=False)
    make_plot(summary)

    print(f"Saved {OUT_TSV}")
    print(f"Saved {OUT_PDF}")
    print(f"Saved {OUT_PNG}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
