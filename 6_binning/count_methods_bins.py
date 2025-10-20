#!/usr/bin/env python3
import pandas as pd
import os

# Hardcoded paths to the .stats files
stats_files = {
    "Integrated": "metawrap_50_10_bins.stats",
    "MetaCC": "binsA.stats",
    "ImputeCC": "binsB.stats",
    "bin3C": "binsC.stats"
}

# Thresholds
completeness_thresholds = [50, 70, 90]
contamination_thresholds = [5, 10]

# Prepare tables
tables = {cont: pd.DataFrame(index=stats_files.keys(), columns=completeness_thresholds) for cont in contamination_thresholds}

for group, filepath in stats_files.items():
    if not os.path.exists(filepath):
        print(f"[WARNING] Missing file: {filepath}")
        continue

    df = pd.read_csv(filepath, sep="\t")

    if 'completeness' not in df.columns or 'contamination' not in df.columns:
        print(f"[ERROR] Missing required columns in: {filepath}")
        continue

    for cont in contamination_thresholds:
        for comp in completeness_thresholds:
            count = df[(df['completeness'] >= comp) & (df['contamination'] <= cont)].shape[0]
            tables[cont].at[group, comp] = count

# Display each table
for cont in contamination_thresholds:
    print(f"\nHuman gut bins with contamination ≤ {cont}")
    print(tables[cont].rename_axis("Group").rename(columns=lambda x: f"Completeness ≥ {x}"))
