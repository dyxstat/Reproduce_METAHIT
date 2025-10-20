#!/usr/bin/env python3
"""
Run Virgo per-contig on QC-passed viral contigs connected to bin.32 and bin.39
Report how many viral contigs are linked (regardless of annotation).
Annotate with Virgo if possible and produce virus_order_contigs.csv with strengths.
"""

import os
import sys
import subprocess
import pandas as pd
from Bio import SeqIO
from scipy.sparse import load_npz

# Paths
OUTDIR = "/work/dulab/Shiyuan/metahit/m4/MGE_output"
GENOMAD_OUT = os.path.join(OUTDIR, "genomad_output")
ASSEMBLY = "combined_contigs"

# Viral contigs (QC-passed from 10_MGE.sh)
VIRAL_FASTA = os.path.join(OUTDIR, "checkv_output", "virus_qc.fna")

# Combined contigs FASTA (for contact matrix index)
COMBINED_FASTA = "/work/dulab/Shiyuan/metahit/m4/reassembly_output/combined/combined_contigs.fa"

# Contact matrix
CONTACT_MATRIX = "/work/dulab/Shiyuan/metahit/m4/contact_output_MGE/denoised_contact_matrix_normcc.npz"

# Virgo setup
VIRGO_REPO = "/work/dulab/Shiyuan/metahit/external/Virgo"
VIRGO_SCRIPT = os.path.join(VIRGO_REPO, "src", "virgo.py")
VIRGO_DB = "/work/dulab/Shiyuan/metahit/external/Virgo_db"
THREADS = "8"

# Safety checks
if not os.path.exists(os.path.join(VIRGO_DB, "database.pkl")):
    sys.exit(f"[ERROR] Virgo DB not found at {VIRGO_DB} (expected database.pkl)")
if not os.path.exists(CONTACT_MATRIX):
    sys.exit(f"[ERROR] Contact matrix not found at {CONTACT_MATRIX}")
if not os.path.exists(COMBINED_FASTA):
    sys.exit(f"[ERROR] Combined FASTA not found at {COMBINED_FASTA}")
if not os.path.exists(VIRAL_FASTA):
    sys.exit(f"[ERROR] Viral QC FASTA not found at {VIRAL_FASTA}")

# Load viral QC contigs
viral_records = {rec.id: rec for rec in SeqIO.parse(VIRAL_FASTA, "fasta")}
viral_ids = set(viral_records.keys())

# Load contact matrix and contig names
mat = load_npz(CONTACT_MATRIX).tocoo()
contigs = [line[1:].split()[0] for line in open(COMBINED_FASTA) if line.startswith(">")]

# Step 1: find viral contigs connected to bin.32 and bin.39
bins_interest = ["bin.32", "bin.39"]  # check prefix only
linked_contigs = {b: set() for b in bins_interest}

for r, c, v in zip(mat.row, mat.col, mat.data):
    if v <= 0:
        continue
    a, b = contigs[r], contigs[c]
    for viral, host in [(a, b), (b, a)]:
        if viral in viral_ids:
            for bprefix in bins_interest:
                if host.startswith(bprefix):
                    linked_contigs[bprefix].add(viral)

print(f"[INFO] Found {len(linked_contigs['bin.32'])} viral contigs linked with bin.32")
print(f"[INFO] Found {len(linked_contigs['bin.39'])} viral contigs linked with bin.39")

# Step 2: write per-contig FASTAs for Virgo
virgo_subset_dir = os.path.join(OUTDIR, "virgo_percontig")
os.makedirs(virgo_subset_dir, exist_ok=True)

percontig_dirs = []
for binname, contigset in linked_contigs.items():
    for contig in contigset:
        subdir = os.path.join(virgo_subset_dir, contig)
        os.makedirs(subdir, exist_ok=True)
        fa_path = os.path.join(subdir, f"{contig}.fasta")
        SeqIO.write(viral_records[contig], fa_path, "fasta")
        percontig_dirs.append((contig, binname, subdir, fa_path))

print(f"[INFO] Prepared {len(percontig_dirs)} viral contigs for Virgo annotation")

# Step 3: run Virgo on each contig
contig_to_order = {}
for contig, binname, subdir, fa_path in percontig_dirs:
    out_ann = os.path.join(subdir, "virgo_annotation")
    os.makedirs(out_ann, exist_ok=True)
    res_file = os.path.join(out_ann, "results.csv")
    if not os.path.exists(res_file):
        print(f"[INFO] Running Virgo on {contig}")
        virgo_cmd = [
            "python", VIRGO_SCRIPT,
            "-d", VIRGO_DB,
            "-i", subdir,
            "-o", out_ann,
            "-t", THREADS
        ]
        subprocess.run(virgo_cmd, check=True)
    if os.path.exists(res_file):
        df = pd.read_csv(res_file)
        if not df.empty and "Order" in df.columns:
            order = str(df.loc[0, "Order"])
            if order.lower() != "nan" and order.strip() != "":
                contig_to_order[contig] = order

print(f"[INFO] Virgo annotated {len(contig_to_order)} viral contigs with order info")

# Step 4: accumulate contact strengths between annotated viral contigs and bins
contig_bin_strength = {}
for r, c, v in zip(mat.row, mat.col, mat.data):
    if v <= 0:
        continue
    a, b = contigs[r], contigs[c]
    for contig, host in [(a, b), (b, a)]:
        if contig in contig_to_order:
            for bprefix in bins_interest:
                if host.startswith(bprefix):
                    bin_id = bprefix  # "bin.32" or "bin.39"
                    key = (contig, bin_id)
                    contig_bin_strength[key] = contig_bin_strength.get(key, 0) + int(v)

# Aggregate by order
order_hits = {}
order_counter = {}
for (contig, bin_id), strength in contig_bin_strength.items():
    order = contig_to_order[contig]
    order_counter.setdefault(order, 0)
    order_counter[order] += 1
    order_tag = f"{order}_{order_counter[order]}"
    order_hits.setdefault(bin_id, {})
    order_hits[bin_id][order_tag] = order_hits[bin_id].get(order_tag, 0) + strength

# Step 5: write output
out_file = os.path.join(virgo_subset_dir, "virus_order_contigs.csv")
with open(out_file, "w") as fh:
    fh.write("bin,contigs\n")
    for binname, d in order_hits.items():
        entries = ",".join(f"{m}({w})" for m, w in d.items())
        fh.write(f"{binname},{entries}\n")

print(f"[INFO] Wrote {out_file}")
