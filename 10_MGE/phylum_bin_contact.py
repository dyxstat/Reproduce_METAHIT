#!/usr/bin/env python3
"""
Count viral and plasmid contigs contacting bins in the top 4 phyla.
Output: phylum_bin_contact.tsv in the same folder as this script.
"""

import os
import pandas as pd
from Bio import SeqIO
from scipy.sparse import load_npz

# -------------------------
# Input files (adjust paths)
# -------------------------
gtdb_summary = "gtdbtk.bac120.summary.tsv"  # put this file in same folder as script
combined_fasta = "/work/dulab/Shiyuan/metahit/m4/reassembly_output/combined/combined_contigs.fa"
contact_matrix = "/work/dulab/Shiyuan/metahit/m4/contact_output_MGE/denoised_contact_matrix_normcc.npz"
viral_fasta = "/work/dulab/Shiyuan/metahit/m4/MGE_output/checkv_output/virus_qc.fna"
plasmid_fasta = "/work/dulab/Shiyuan/metahit/m4/MGE_output/genomad_output/combined_contigs_summary/combined_contigs_plasmid.fna"

# Top 4 phyla of interest
target_phyla = {"Bacteroidota","Bacillota_A","Pseudomonadota","Actinomycetota"}

# -------------------------
# Load GTDB bin -> phylum map
# -------------------------
gtdb = pd.read_csv(gtdb_summary, sep="\t")
bin_to_phylum = {}
for _, row in gtdb.iterrows():
    bin_name = row["user_genome"].replace(".fa","")
    lineage = row["classification"].split(";")
    phylum = [x[3:] for x in lineage if x.startswith("p__")]
    if phylum:
        bin_to_phylum[bin_name] = phylum[0]

# -------------------------
# Load contigs + sets
# -------------------------
contigs = [line[1:].split()[0] for line in open(combined_fasta) if line.startswith(">")]
viral_ids = {rec.id for rec in SeqIO.parse(viral_fasta,"fasta")}
plasmid_ids = {rec.id for rec in SeqIO.parse(plasmid_fasta,"fasta")}
mat = load_npz(contact_matrix).tocoo()

# -------------------------
# Count contacts per phylum
# -------------------------
stats = {p: {"Viral":0,"Plasmid":0,"Bins":set()} for p in target_phyla}

for r,c,v in zip(mat.row, mat.col, mat.data):
    if v <= 0: continue
    a,b = contigs[r], contigs[c]

    for query, host in [(a,b),(b,a)]:
        if host.startswith("bin"):
            bin_id = host.split("_",1)[0]  # e.g. bin.32
            phylum = bin_to_phylum.get(bin_id)
            if phylum in target_phyla:
                stats[phylum]["Bins"].add(bin_id)
                if query in viral_ids:
                    stats[phylum]["Viral"] += 1
                elif query in plasmid_ids:
                    stats[phylum]["Plasmid"] += 1

# -------------------------
# Write output table
# -------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
out_file = os.path.join(script_dir,"phylum_bin_contact.tsv")

with open(out_file,"w") as fh:
    fh.write("Phylum\tViral\tPlasmid\tBins\n")
    for p in target_phyla:
        v = stats[p]["Viral"]
        pl = stats[p]["Plasmid"]
        bins = ",".join(sorted(stats[p]["Bins"])) if stats[p]["Bins"] else "NA"
        fh.write(f"{p}\t{v}\t{pl}\t{bins}\n")

print(f"[INFO] Wrote results to {out_file}")
