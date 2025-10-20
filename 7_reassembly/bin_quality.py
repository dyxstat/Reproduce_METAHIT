#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# ---------------------------------------------------------------------
# Font settings for main figure
# ---------------------------------------------------------------------
mpl.rcParams.update({
    'font.size': 30,
    'axes.titlesize': 50,
    'axes.labelsize': 30,
    'legend.fontsize': 30,
    'xtick.labelsize': 30,
    'ytick.labelsize': 30,
})

# ---------------------------------------------------------------------
# File paths
# ---------------------------------------------------------------------
REFINE_STATS = "metawrap_50_10_bins.stats"
REASSEMBLY_STATS = "reassembled_bins.stats"

# Load and filter
refine_df = pd.read_csv(REFINE_STATS, sep="\t")
reassembly_df = pd.read_csv(REASSEMBLY_STATS, sep="\t")

# Filter bins: ≥50% completeness and ≤10% contamination
refine = refine_df[(refine_df['completeness'] >= 50) & (refine_df['contamination'] <= 10)].copy()
reassembly = reassembly_df[(reassembly_df['completeness'] >= 50) & (reassembly_df['contamination'] <= 10)].copy()

# ------- Top plot: sort by completeness -------
refine_comp_sorted = refine.sort_values(by='completeness', ascending=False).reset_index(drop=True)
reassembly_comp_sorted = reassembly.sort_values(by='completeness', ascending=False).reset_index(drop=True)
min_len_comp = min(len(refine_comp_sorted), len(reassembly_comp_sorted))
refine_comp_sorted = refine_comp_sorted.iloc[:min_len_comp]
reassembly_comp_sorted = reassembly_comp_sorted.iloc[:min_len_comp]

# ------- Bottom plot: sort by contamination -------
refine_cont_sorted = refine.sort_values(by='contamination').reset_index(drop=True)
reassembly_cont_sorted = reassembly.sort_values(by='contamination').reset_index(drop=True)
min_len_cont = min(len(refine_cont_sorted), len(reassembly_cont_sorted))
refine_cont_sorted = refine_cont_sorted.iloc[:min_len_cont]
reassembly_cont_sorted = reassembly_cont_sorted.iloc[:min_len_cont]

# ------- Plotting main figure -------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=False)

before_color = '#1f77b4'  # blue
after_color = '#ff7f0e'   # orange

# Completion curve
ax1.plot(refine_comp_sorted['completeness'].values, label='Before reassembly', color=before_color)
ax1.plot(reassembly_comp_sorted['completeness'].values, label='After reassembly', color=after_color)
ax1.set_ylabel("Completion")
ax1.set_title("Human")
ax1.set_ylim(50, 100)

# Contamination curve
ax2.plot(refine_cont_sorted['contamination'].values, label='Before reassembly', color=before_color)
ax2.plot(reassembly_cont_sorted['contamination'].values, label='After reassembly', color=after_color)
ax2.set_ylabel("Contamination")
ax2.set_ylim(0, 10)

fig.align_ylabels([ax1, ax2])
ax1.set_position([0.15, 0.55, 0.8, 0.35])
ax2.set_position([0.15, 0.1, 0.8, 0.35])

plt.savefig("bin_quality_human.pdf")
plt.show()  # <--- show both plots instead of closing

# ---------------------------------------------------------------------
# Standalone legend (block style, same as fig3_legend.pdf)
# ---------------------------------------------------------------------
def create_bin_quality_legend(output_file="bin_quality_legend.pdf"):
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': 30,
        'legend.fontsize': 30
    })

    fig, ax = plt.subplots(figsize=(4, 1))
    ax.axis('off')

    before_color = '#1f77b4'  # blue
    after_color = '#ff7f0e'   # orange

    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, color=before_color, alpha=0.8, label='Before reassembly'),
        plt.Rectangle((0, 0), 1, 1, color=after_color, alpha=0.8, label='After reassembly')
    ]

    legend = fig.legend(handles=legend_elements, loc='center', ncol=2, 
                        fontsize=30, frameon=False)
    legend.get_frame().set_edgecolor('none')
    legend.get_frame().set_linewidth(0)

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Legend saved as '{output_file}'")
    plt.close('all')

# Create the standalone legend
create_bin_quality_legend()
