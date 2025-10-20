#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reassembly Metrics Evaluation Script

Calculates 7 key metrics to evaluate reassembly effectiveness using population-level averages:
1. Average Completeness (CheckM2)
2. Average Contamination (CheckM2) 
3. Average N50 (assembly contiguity)
4. Average L50 (number of contigs needed for N50)
5. Average Total Assembly Size (total bp per bin)
6. Average Number of Contigs (fragmentation measure)
7. Average Largest Contig Size (longest single sequence)

Author: Metahit Pipeline
"""

import argparse
import os
import sys
import subprocess
import pandas as pd
import numpy as np
import glob
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate reassembly effectiveness - 7 key metrics using population-level averages"
    )
    parser.add_argument('--original_bins', required=True, 
                       help='Directory containing original bins (*.fa)')
    parser.add_argument('--reassembled_bins', required=True,
                       help='Directory containing reassembled bins (*.fa)')
    parser.add_argument('--outdir', required=True,
                       help='Output directory for evaluation results')
    parser.add_argument('-t', '--threads', default='8',
                       help='Number of threads (default: 8)')
    parser.add_argument('--skip_checkm2', action='store_true',
                       help='Skip CheckM2 analysis (assembly stats only)')
    return parser.parse_args()


def run_command(cmd, check=True):
    """Run shell command with error handling"""
    print(f"[CMD] {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, check=check, 
                              capture_output=True, text=True)
        if result.stdout.strip():
            print(f"[STDOUT] {result.stdout.strip()}")
        return result
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {cmd}")
        if e.stderr:
            print(f"[STDERR] {e.stderr}")
        if check:
            raise
        return e


def calculate_assembly_stats(fasta_file):
    """Calculate N50, L50, and other assembly statistics"""
    if not os.path.exists(fasta_file):
        return None
    
    lengths = []
    total_length = 0
    
    # Read FASTA and get contig lengths
    with open(fasta_file, 'r') as f:
        current_length = 0
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_length > 0:
                    lengths.append(current_length)
                    total_length += current_length
                current_length = 0
            else:
                current_length += len(line)
        
        # Don't forget the last contig
        if current_length > 0:
            lengths.append(current_length)
            total_length += current_length
    
    if not lengths:
        return None
    
    # Sort lengths in descending order
    lengths.sort(reverse=True)
    
    # Calculate N50 and L50
    half_total = total_length / 2
    cumulative = 0
    n50 = 0
    l50 = 0
    
    for i, length in enumerate(lengths):
        cumulative += length
        if cumulative >= half_total:
            n50 = length
            l50 = i + 1
            break
    
    stats = {
        'total_length': total_length,
        'num_contigs': len(lengths),
        'n50': n50,
        'l50': l50,
        'largest_contig': max(lengths),
        'mean_length': np.mean(lengths),
        'median_length': np.median(lengths)
    }
    
    return stats


def run_checkm2(bin_dir, output_dir, threads):
    """Run CheckM2 analysis following metaWRAP approach exactly"""
    os.makedirs(output_dir, exist_ok=True)
    
    print("[INFO] Running CheckM2...")
    checkm_output = os.path.join(output_dir, 'checkm2_output')
    if os.path.exists(checkm_output):
        run_command(f"rm -r {checkm_output}")
    os.makedirs(checkm_output, exist_ok=True)
    
    # Create temporary directory for CheckM2
    tmp_dir = os.path.join(output_dir, 'checkm2_tmp')
    os.makedirs(tmp_dir, exist_ok=True)
    
    # Run CheckM2 directly (assuming checkm2 environment is active, like metaWRAP does)
    result = run_command(f"checkm2 predict -i {bin_dir} -o {checkm_output} -x fa -t {threads} --tmpdir {tmp_dir}", check=False)
    
    # Clean up temp directory
    if os.path.exists(tmp_dir):
        run_command(f"rm -r {tmp_dir}")
    
    if result.returncode != 0:
        print("[WARNING] CheckM2 analysis failed")
        print(f"[ERROR] Return code: {result.returncode}")
        if result.stderr:
            print(f"[ERROR] STDERR: {result.stderr}")
        return None
    
    # Parse CheckM2 results
    results_file = os.path.join(checkm_output, 'quality_report.tsv')
    if os.path.exists(results_file):
        df = pd.read_csv(results_file, sep='\t')
        return df
    
    return None


def analyze_bins(bin_dir, output_dir, name_prefix, threads, skip_checkm2=False):
    """Analyze a set of bins and return population-level statistics"""
    print(f"[INFO] Analyzing bins in {bin_dir}...")
    
    # Create output subdirectory
    analysis_dir = os.path.join(output_dir, f"{name_prefix}_analysis")
    os.makedirs(analysis_dir, exist_ok=True)
    
    # Find all FASTA files (metaWRAP uses .fa extension)
    fasta_files = glob.glob(os.path.join(bin_dir, "*.fa"))
    
    if not fasta_files:
        print(f"[WARNING] No .fa files found in {bin_dir}")
        return {}
    
    print(f"[INFO] Found {len(fasta_files)} bins to analyze")
    
    # Calculate assembly statistics for each bin
    all_stats = []
    for fasta_file in fasta_files:
        bin_name = os.path.splitext(os.path.basename(fasta_file))[0]
        print(f"[INFO] Processing {bin_name}...")
        
        stats = calculate_assembly_stats(fasta_file)
        if stats is not None:
            all_stats.append(stats)
    
    if not all_stats:
        print(f"[WARNING] No valid bins found in {bin_dir}")
        return {}
    
    # Calculate population-level averages for assembly metrics
    population_stats = {
        'num_bins': len(all_stats),
        'avg_total_length': np.mean([s['total_length'] for s in all_stats]),
        'avg_num_contigs': np.mean([s['num_contigs'] for s in all_stats]),
        'avg_n50': np.mean([s['n50'] for s in all_stats]),
        'avg_l50': np.mean([s['l50'] for s in all_stats]),
        'avg_largest_contig': np.mean([s['largest_contig'] for s in all_stats]),
        'avg_mean_length': np.mean([s['mean_length'] for s in all_stats]),
        'avg_median_length': np.mean([s['median_length'] for s in all_stats])
    }
    
    # Run CheckM2 analysis if not skipped - activate checkm2 environment first
    if not skip_checkm2 and len(fasta_files) > 0:
        print(f"[INFO] Switching to CheckM2 environment...")
        
        # Switch to checkm2 environment like metaWRAP does
        run_command("conda deactivate")
        run_command("conda activate checkm2")
        
        print(f"[INFO] Running CheckM2 analysis on {len(fasta_files)} bins...")
        checkm_df = run_checkm2(bin_dir, analysis_dir, threads)
        
        # Switch back to metahit_env like metaWRAP does
        run_command("conda deactivate")
        run_command("conda activate metahit_env")
        
        if checkm_df is not None:
            # Calculate population-level averages for CheckM2 metrics
            population_stats['avg_completeness'] = checkm_df['Completeness'].mean()
            population_stats['avg_contamination'] = checkm_df['Contamination'].mean()
            print(f"[INFO] CheckM2 results: Avg completeness = {population_stats['avg_completeness']:.2f}%, Avg contamination = {population_stats['avg_contamination']:.2f}%")
        else:
            print("[WARNING] CheckM2 analysis failed - continuing with assembly statistics only")
    else:
        print("[INFO] CheckM2 analysis skipped")
    
    return population_stats


def compare_population_metrics(original_stats, reassembled_stats, output_dir):
    """Compare population-level metrics between original and reassembled bins"""
    print("[INFO] Comparing population-level metrics...")
    
    if not original_stats or not reassembled_stats:
        print("[ERROR] Cannot compare - missing statistics for one or both bin sets")
        return {}
    
    comparison = {
        'original_num_bins': original_stats['num_bins'],
        'reassembled_num_bins': reassembled_stats['num_bins'],
    }
    
    # Compare the 7 key metrics
    metrics = [
        ('completeness', 'avg_completeness'),
        ('contamination', 'avg_contamination'), 
        ('n50', 'avg_n50'),
        ('l50', 'avg_l50'),
        ('total_length', 'avg_total_length'),
        ('num_contigs', 'avg_num_contigs'),
        ('largest_contig', 'avg_largest_contig')
    ]
    
    for metric_name, stat_key in metrics:
        if stat_key in original_stats and stat_key in reassembled_stats:
            original_val = original_stats[stat_key]
            reassembled_val = reassembled_stats[stat_key]
            
            comparison[f'original_{metric_name}'] = original_val
            comparison[f'reassembled_{metric_name}'] = reassembled_val
            comparison[f'{metric_name}_change'] = reassembled_val - original_val
            
            if metric_name in ['n50', 'largest_contig']:
                # Calculate fold change for metrics where bigger is better
                comparison[f'{metric_name}_fold_change'] = reassembled_val / original_val if original_val > 0 else float('inf')
            elif metric_name == 'num_contigs':
                # For contig count, negative change is improvement (fewer contigs)
                comparison[f'{metric_name}_improvement'] = original_val - reassembled_val
    
    # Save detailed comparison
    comparison_file = os.path.join(output_dir, 'population_comparison.tsv')
    comparison_df = pd.DataFrame([comparison])
    comparison_df.to_csv(comparison_file, sep='\t', index=False)
    print(f"[INFO] Population comparison saved to {comparison_file}")
    
    return comparison


def generate_summary(comparison, original_stats, reassembled_stats, output_dir):
    """Generate summary report for population-level comparison"""
    print("[INFO] Generating summary report...")
    
    summary_file = os.path.join(output_dir, 'reassembly_population_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("REASSEMBLY EVALUATION SUMMARY - POPULATION-LEVEL METRICS\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("OVERVIEW:\n")
        f.write(f"Original bins analyzed: {comparison['original_num_bins']}\n")
        f.write(f"Reassembled bins analyzed: {comparison['reassembled_num_bins']}\n\n")
        
        # 1. Completeness
        f.write("1. AVERAGE COMPLETENESS:\n")
        if 'original_completeness' in comparison:
            f.write(f"   Original: {comparison['original_completeness']:.2f}%\n")
            f.write(f"   Reassembled: {comparison['reassembled_completeness']:.2f}%\n")
            f.write(f"   Change: {comparison['completeness_change']:+.2f}%\n")
            improvement = "IMPROVED" if comparison['completeness_change'] > 0 else "DECREASED" if comparison['completeness_change'] < 0 else "UNCHANGED"
            f.write(f"   Result: {improvement}\n\n")
        else:
            f.write("   No completeness data available (CheckM2 skipped)\n\n")
        
        # 2. Contamination  
        f.write("2. AVERAGE CONTAMINATION:\n")
        if 'original_contamination' in comparison:
            f.write(f"   Original: {comparison['original_contamination']:.2f}%\n")
            f.write(f"   Reassembled: {comparison['reassembled_contamination']:.2f}%\n")
            f.write(f"   Change: {comparison['contamination_change']:+.2f}%\n")
            improvement = "IMPROVED" if comparison['contamination_change'] < 0 else "WORSENED" if comparison['contamination_change'] > 0 else "UNCHANGED"
            f.write(f"   Result: {improvement}\n\n")
        else:
            f.write("   No contamination data available (CheckM2 skipped)\n\n")
        
        # 3. N50
        f.write("3. AVERAGE N50 (Assembly Contiguity):\n")
        f.write(f"   Original: {comparison['original_n50']:.0f} bp\n")
        f.write(f"   Reassembled: {comparison['reassembled_n50']:.0f} bp\n")
        f.write(f"   Change: {comparison['n50_change']:+.0f} bp\n")
        f.write(f"   Fold change: {comparison['n50_fold_change']:.2f}x\n")
        improvement = "IMPROVED" if comparison['n50_fold_change'] > 1 else "DECREASED" if comparison['n50_fold_change'] < 1 else "UNCHANGED"
        f.write(f"   Result: {improvement}\n\n")
        
        # 4. L50
        f.write("4. AVERAGE L50 (Contigs for N50):\n")
        f.write(f"   Original: {comparison['original_l50']:.1f}\n")
        f.write(f"   Reassembled: {comparison['reassembled_l50']:.1f}\n")
        f.write(f"   Change: {comparison['l50_change']:+.1f}\n")
        improvement = "IMPROVED" if comparison['l50_change'] < 0 else "WORSENED" if comparison['l50_change'] > 0 else "UNCHANGED"
        f.write(f"   Result: {improvement}\n\n")
        
        # 5. Total assembly size
        f.write("5. AVERAGE TOTAL ASSEMBLY SIZE:\n")
        f.write(f"   Original: {comparison['original_total_length']:.0f} bp\n")
        f.write(f"   Reassembled: {comparison['reassembled_total_length']:.0f} bp\n")
        f.write(f"   Change: {comparison['total_length_change']:+.0f} bp\n")
        improvement = "INCREASED" if comparison['total_length_change'] > 0 else "DECREASED" if comparison['total_length_change'] < 0 else "UNCHANGED"
        f.write(f"   Result: {improvement}\n\n")
        
        # 6. Number of contigs
        f.write("6. AVERAGE NUMBER OF CONTIGS:\n")
        f.write(f"   Original: {comparison['original_num_contigs']:.1f}\n")
        f.write(f"   Reassembled: {comparison['reassembled_num_contigs']:.1f}\n")
        f.write(f"   Change: {comparison['num_contigs_change']:+.1f}\n")
        f.write(f"   Improvement: {comparison['num_contigs_improvement']:+.1f} (fewer is better)\n")
        improvement = "IMPROVED" if comparison['num_contigs_improvement'] > 0 else "WORSENED" if comparison['num_contigs_improvement'] < 0 else "UNCHANGED"
        f.write(f"   Result: {improvement}\n\n")
        
        # 7. Largest contig
        f.write("7. AVERAGE LARGEST CONTIG SIZE:\n")
        f.write(f"   Original: {comparison['original_largest_contig']:.0f} bp\n")
        f.write(f"   Reassembled: {comparison['reassembled_largest_contig']:.0f} bp\n")
        f.write(f"   Change: {comparison['largest_contig_change']:+.0f} bp\n")
        f.write(f"   Fold change: {comparison['largest_contig_fold_change']:.2f}x\n")
        improvement = "IMPROVED" if comparison['largest_contig_fold_change'] > 1 else "DECREASED" if comparison['largest_contig_fold_change'] < 1 else "UNCHANGED"
        f.write(f"   Result: {improvement}\n\n")
        
        # Overall assessment
        f.write("OVERALL ASSESSMENT:\n")
        improvements = 0
        total_metrics = 7
        
        if 'completeness_change' in comparison and comparison['completeness_change'] > 0:
            improvements += 1
        if 'contamination_change' in comparison and comparison['contamination_change'] < 0:
            improvements += 1
        if comparison['n50_fold_change'] > 1:
            improvements += 1
        if comparison['l50_change'] < 0:
            improvements += 1
        if comparison['total_length_change'] > 0:
            improvements += 1
        if comparison['num_contigs_improvement'] > 0:
            improvements += 1
        if comparison['largest_contig_fold_change'] > 1:
            improvements += 1
            
        if 'completeness_change' not in comparison:
            total_metrics = 5  # No CheckM2 data
            
        f.write(f"Metrics improved: {improvements}/{total_metrics}\n")
        
        if improvements >= total_metrics * 0.7:
            f.write("CONCLUSION: Reassembly was HIGHLY SUCCESSFUL\n")
        elif improvements >= total_metrics * 0.5:
            f.write("CONCLUSION: Reassembly was MODERATELY SUCCESSFUL\n")
        elif improvements > 0:
            f.write("CONCLUSION: Reassembly showed SOME IMPROVEMENT\n")
        else:
            f.write("CONCLUSION: Reassembly showed NO IMPROVEMENT\n")
    
    print(f"[INFO] Summary saved to {summary_file}")
    return summary_file


def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    print(f"[INFO] Starting reassembly evaluation (population-level metrics)...")
    print(f"[INFO] Original bins: {args.original_bins}")
    print(f"[INFO] Reassembled bins: {args.reassembled_bins}")
    print(f"[INFO] Output directory: {args.outdir}")
    print(f"[INFO] Threads: {args.threads}")
    print(f"[INFO] Skip CheckM2: {args.skip_checkm2}")
    
    # Analyze original bins
    print("\n[INFO] Analyzing original bins...")
    original_stats = analyze_bins(
        args.original_bins, args.outdir, "original", args.threads, args.skip_checkm2
    )
    
    if not original_stats:
        print("[ERROR] No original bins found or analyzed. Check the original bins directory.")
        return
    
    print(f"[INFO] Original bins summary: {original_stats['num_bins']} bins analyzed")
    
    # Analyze reassembled bins
    print("\n[INFO] Analyzing reassembled bins...")
    reassembled_stats = analyze_bins(
        args.reassembled_bins, args.outdir, "reassembled", args.threads, args.skip_checkm2
    )
    
    if not reassembled_stats:
        print("[ERROR] No reassembled bins found or analyzed. Check the reassembled bins directory.")
        return
    
    print(f"[INFO] Reassembled bins summary: {reassembled_stats['num_bins']} bins analyzed")
    
    # Compare population-level metrics
    print(f"\n[INFO] Comparing population-level metrics...")
    comparison = compare_population_metrics(original_stats, reassembled_stats, args.outdir)
    
    if not comparison:
        print("[ERROR] Failed to generate comparison")
        return
    
    # Generate summary report
    summary_file = generate_summary(comparison, original_stats, reassembled_stats, args.outdir)
    
    print("\n" + "=" * 70)
    print("REASSEMBLY EVALUATION COMPLETED - POPULATION-LEVEL METRICS")
    print("=" * 70)
    print(f"Original bins: {comparison['original_num_bins']}")
    print(f"Reassembled bins: {comparison['reassembled_num_bins']}")
    
    if 'completeness_change' in comparison:
        print(f"Completeness change: {comparison['completeness_change']:+.2f}%")
        print(f"Contamination change: {comparison['contamination_change']:+.2f}%")
    else:
        print("CheckM2 analysis: SKIPPED")
    
    print(f"N50 fold change: {comparison['n50_fold_change']:.2f}x")
    print(f"Contig count improvement: {comparison['num_contigs_improvement']:+.1f}")
    print(f"Largest contig fold change: {comparison['largest_contig_fold_change']:.2f}x")
    
    print(f"\nDetailed results saved to: {args.outdir}")
    print("Key files:")
    print(f"  - Population comparison: {args.outdir}/population_comparison.tsv")
    print(f"  - Summary report: {summary_file}")


if __name__ == '__main__':
    main()