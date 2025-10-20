#!/usr/bin/env python3
import pandas as pd

# File paths
REFINE_STATS = "metawrap_50_10_bins.stats"
REASSEMBLY_STATS = "reassembled_bins.stats"

def calculate_avg_completeness(stats_file):
    """Calculate average completeness for high-quality bins."""
    df = pd.read_csv(stats_file, sep="\t")
    
    # Filter bins: ≥50% completeness and ≤10% contamination
    filtered_df = df[(df['completeness'] >= 50) & (df['contamination'] <= 10)]
    
    # Calculate average completeness
    avg_comp = filtered_df['completeness'].mean()
    
    return avg_comp

def main():
    # Calculate for both datasets
    refine_avg_comp = calculate_avg_completeness(REFINE_STATS)
    reassembly_avg_comp = calculate_avg_completeness(REASSEMBLY_STATS)
    
    print(f"Average completeness before reassembly: {refine_avg_comp:.2f}%")
    print(f"Average completeness after reassembly: {reassembly_avg_comp:.2f}%")

if __name__ == "__main__":
    main()