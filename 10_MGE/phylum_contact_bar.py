#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

# Input/output
file_in = "phylum_bin_contact.tsv"
file_out = "phylum_contact_bar.pdf"

# Load data
df = pd.read_csv(file_in, sep="\t")

# Compute total and sort descending
df["Total"] = df["Viral"] + df["Plasmid"]
df = df.sort_values("Total", ascending=False)

# Data
phylums = df["Phylum"].tolist()
viral = df["Viral"].astype(int).tolist()
plasmid = df["Plasmid"].astype(int).tolist()

# Plot
fig, ax = plt.subplots(figsize=(20, 12))  # enlarged from (12, 8)
bar_width = 0.6
x = range(len(phylums))

# Bars
p1 = ax.bar(x, plasmid, color="#ffcccc", width=bar_width, label="Plasmid contigs")
p2 = ax.bar(x, viral, bottom=plasmid, color="#cce6ff", width=bar_width, label="Viral contigs")

# X-axis formatting
ax.set_xticks(x)
ax.set_xticklabels(phylums, fontsize=36, ha="center")  # doubled from 20

# Y-axis label
ax.set_ylabel("Number of contigs", fontsize=40)  # doubled from 20

# Tick font sizes
ax.tick_params(axis='y', labelsize=36)  # doubled from ~18

# Legend at bottom
fig.legend(
    handles=[p1, p2],   # Keep order: Plasmid then Viral
    loc="lower center",
    ncol=2,
    fontsize=36,  # doubled from 18
    frameon=False
)

# Adjust layout to reduce space for legend
plt.tight_layout(rect=[0, 0.07, 1, 1])

# Save
plt.savefig(file_out, dpi=300)
plt.close()

print(f"Saved stacked barplot to {file_out}")
