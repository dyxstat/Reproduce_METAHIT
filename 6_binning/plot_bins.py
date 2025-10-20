#!/usr/bin/env python3
import pandas as pd
import os
import matplotlib.pyplot as plt

# Mapping of binning methods to stats files
stats_files = {
    "METAHIT": "metawrap_50_10_bins.stats",
    "bin3C": "binsC.stats",
    "MetaCC": "binsA.stats",
    "ImputeCC": "binsB.stats"
}

# Contamination thresholds and output filenames
contamination_thresholds = {
    5: "bin_counts_human_5.pdf",
    10: "bin_counts_human_10.pdf"
}

# Colors: light to dark (for plotting)
plot_colors = ['#a3b18a', '#588157', '#3a5a40']
plot_thresholds = [50, 70, 90]
plot_labels = ["Completeness ≥ 50", "Completeness ≥ 70", "Completeness ≥ 90"]

# Colors and labels reversed for legend (dark to light)
legend_colors = plot_colors[::-1]
legend_labels = plot_labels[::-1]

# Set global font sizes via rcParams
plt.rcParams.update({
    'font.size': 30,             # Base font size (for ticks)
    'axes.labelsize': 30,        # X and Y label size
    'axes.titlesize': 50,        # Title font size
    'xtick.labelsize': 30,       # Tick label size
    'ytick.labelsize': 30,
    'legend.fontsize': 30        # Legend font size
})

# Process each contamination threshold
for cont, output_file in contamination_thresholds.items():
    table = pd.DataFrame(index=stats_files.keys(), columns=plot_labels)

    for group, filepath in stats_files.items():
        if not os.path.exists(filepath):
            print(f"[WARNING] Missing file: {filepath}")
            continue

        df = pd.read_csv(filepath, sep="\t")
        if 'completeness' not in df.columns or 'contamination' not in df.columns:
            print(f"[ERROR] Missing required columns in: {filepath}")
            continue

        filtered = df[df['contamination'] <= cont]
        for comp, label in zip(plot_thresholds, plot_labels):
            count = filtered[filtered['completeness'] >= comp].shape[0]
            table.at[group, label] = count

    df_plot = table.fillna(0).astype(int)
    df_plot = df_plot.loc[["ImputeCC", "MetaCC", "bin3C", "METAHIT"]]  # reverse for top-down

    fig, ax = plt.subplots(figsize=(12, 7))  # Increased size for better spacing

    # Draw bars in light → dark (so darker overlays shorter lighter ones)
    for color, label in zip(plot_colors, plot_labels):
        ax.barh(df_plot.index, df_plot[label], color=color, label=label)
        
    ax.set_position([0.2, 0.15, 0.7, 0.7])

    #ax.set_xlabel("Number of bins")
    #ax.set_ylabel("Binning Method")
    ax.set_title(f"Human Gut Dataset", weight='bold')

    # # Legend (dark → light for natural semantic reading)
    # handles, _ = ax.get_legend_handles_labels()
    # ax.legend(
    #     handles[::-1], legend_labels,
    #     loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False
    # )

    fig.savefig(output_file, bbox_inches='tight')
    plt.show()
