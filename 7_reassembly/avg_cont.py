#!/usr/bin/env python3
import pandas as pd

# File paths
REFINE_STATS = "metawrap_50_10_bins.stats"
REASSEMBLY_STATS = "reassembled_bins.stats"

def calculate_avg_contamination(stats_file, label):
    """Calculate average contamination for high-quality bins."""
    df = pd.read_csv(stats_file, sep="\t")
    
    # Filter bins: ≥50% completeness and ≤10% contamination
    filtered_df = df[(df['completeness'] >= 50) & (df['contamination'] <= 10)]
    
    # Calculate average contamination
    avg_cont = filtered_df['contamination'].mean()
    
    print(f"{label}:")
    print(f"  Total high-quality bins: {len(filtered_df)}")
    print(f"  Average contamination: {avg_cont:.2f}%")
    
    return avg_cont, len(filtered_df)

def main():
    print("Calculating average contamination scores for high-quality bins")
    print("(≥50% completeness and ≤10% contamination)")
    print("-" * 60)
    
    # Calculate for both datasets
    refine_avg, refine_count = calculate_avg_contamination(REFINE_STATS, "Refined bins")
    print()
    reassembly_avg, reassembly_count = calculate_avg_contamination(REASSEMBLY_STATS, "Reassembled bins")
    
    print()
    print("-" * 60)
    print("SUMMARY:")
    print(f"Decrease in average contamination: {refine_avg - reassembly_avg:.2f} percentage points")
    print(f"Relative improvement: {((refine_avg - reassembly_avg) / refine_avg * 100):.1f}%")

if __name__ == "__main__":
    main()