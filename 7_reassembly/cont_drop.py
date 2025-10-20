#!/usr/bin/env python3
"""
Compare contamination levels between original and reassembled bins
"""

import pandas as pd
import re
import os

def extract_bin_number(bin_name):
    """Extract the bin number from bin name"""
    match = re.match(r'bin\.(\d+)', bin_name)
    return int(match.group(1)) if match else None

def compare_contamination(original_stats, reassembled_stats):
    """Compare contamination between original and reassembled bins"""
    
    # Load the stats files
    original_df = pd.read_csv(original_stats, sep='\t')
    reassembled_df = pd.read_csv(reassembled_stats, sep='\t')
    
    # Extract bin numbers
    original_df['bin_num'] = original_df['bin'].apply(extract_bin_number)
    reassembled_df['bin_num'] = reassembled_df['bin'].apply(extract_bin_number)
    
    # Remove any rows where bin number couldn't be extracted
    original_df = original_df[original_df['bin_num'].notna()]
    reassembled_df = reassembled_df[reassembled_df['bin_num'].notna()]
    
    # Find matching bins (bins that exist in both datasets)
    original_bins = set(original_df['bin_num'])
    reassembled_bins = set(reassembled_df['bin_num'])
    common_bins = original_bins & reassembled_bins
    
    print(f"Original bins: {len(original_bins)}")
    print(f"Reassembled bins: {len(reassembled_bins)}")
    print(f"Common bins: {len(common_bins)}")
    
    if len(common_bins) == 0:
        print("No matching bins found!")
        return 0, 0, 0.0
    
    # Filter to common bins only
    orig_filtered = original_df[original_df['bin_num'].isin(common_bins)].copy()
    reass_filtered = reassembled_df[reassembled_df['bin_num'].isin(common_bins)].copy()
    
    # Sort by bin number for proper matching
    orig_filtered = orig_filtered.sort_values('bin_num').reset_index(drop=True)
    reass_filtered = reass_filtered.sort_values('bin_num').reset_index(drop=True)
    
    # Count bins with lower contamination
    lower_contamination_count = 0
    
    print("\nBin-by-bin comparison:")
    print("Bin\tOriginal\tReassembled\tLower?")
    
    for i, bin_num in enumerate(sorted(common_bins)):
        orig_cont = orig_filtered[orig_filtered['bin_num'] == bin_num]['contamination'].iloc[0]
        reass_cont = reass_filtered[reass_filtered['bin_num'] == bin_num]['contamination'].iloc[0]
        
        is_lower = reass_cont < orig_cont
        if is_lower:
            lower_contamination_count += 1
        
        print(f"{bin_num}\t{orig_cont:.2f}\t\t{reass_cont:.2f}\t\t{'Yes' if is_lower else 'No'}")
    
    total_bins = len(common_bins)
    percentage = (lower_contamination_count / total_bins) * 100
    
    return lower_contamination_count, total_bins, percentage

def main():
    # Check if we're in a stats directory
    if not (os.path.exists('metawrap_50_10_bins.stats') and os.path.exists('reassembled_bins.stats')):
        print("Error: This script should be run in a directory containing both:")
        print("  - metawrap_50_10_bins.stats")
        print("  - reassembled_bins.stats")
        return
    
    print("Comparing contamination levels...")
    print("="*50)
    
    lower_count, total_count, percentage = compare_contamination(
        'metawrap_50_10_bins.stats', 
        'reassembled_bins.stats'
    )
    
    print("\n" + "="*50)
    print("RESULTS:")
    print(f"Bins with lower contamination: {lower_count}")
    print(f"Total comparable bins: {total_count}")
    print(f"Percentage: {percentage:.1f}%")
    print("="*50)

if __name__ == "__main__":
    main()