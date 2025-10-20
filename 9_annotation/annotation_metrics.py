#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotation Metrics Script

Generates a comprehensive table combining:
1. Bin name
2. Bin coverage (weighted average from contig coverage)
3. Completeness (CheckM2)
4. Contamination (CheckM2)
5. Taxonomic classification split into levels (d, p, c, o, f, g, s)
"""

import argparse
import os
import sys
import subprocess
import pandas as pd
import numpy as np
import glob
from pathlib import Path
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate annotation metrics table combining coverage, quality, and taxonomy"
    )
    parser.add_argument('--gtdbtk_dir', required=True,
                       help='Directory containing GTDB-Tk output files')
    parser.add_argument('--bins_dir', required=True,
                       help='Directory containing bin files (*.fa)')
    parser.add_argument('--coverage_file', required=True,
                       help='Coverage file from jgi_summarize_bam_contig_depths')
    parser.add_argument('--outdir', required=True,
                       help='Output directory for results')
    parser.add_argument('--script_dir', required=True,
                       help='Directory containing helper scripts (run_checkm2.sh)')
    parser.add_argument('--threads', default='8',
                       help='Number of threads for CheckM2 (default: 8)')
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


def parse_coverage_file(coverage_file):
    """Parse coverage.txt file to get contig coverage data"""
    print(f"[INFO] Parsing coverage file: {coverage_file}")
    
    coverage_data = {}
    
    with open(coverage_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            
            # Extract contig name (remove the extra info after contig name)
            contig_name_full = parts[0]
            # Try multiple patterns to extract contig name
            
            # Pattern 1: "k141_179321 flag=1 multi=4.0000 len=1233" -> "k141_179321"
            # Pattern 2: "NODE_106_length_522_cov_3.548180" -> "NODE_106_length_522_cov_3.548180"
            
            # First try splitting by space (for k141 pattern)
            contig_name = contig_name_full.split()[0]
            
            # Also store the full name in case bins use full names
            if contig_name_full != contig_name:
                # Store both versions to handle different naming conventions
                contig_names_to_try = [contig_name, contig_name_full]
            else:
                contig_names_to_try = [contig_name]
            
            try:
                contig_len = int(float(parts[1]))
                total_avg_depth = float(parts[2])
                
                # Store coverage data for all possible contig name variants
                for name in contig_names_to_try:
                    coverage_data[name] = {
                        'length': contig_len,
                        'coverage': total_avg_depth
                    }
                    
            except (ValueError, IndexError):
                print(f"[WARNING] Could not parse line: {line}")
                continue
    
    print(f"[INFO] Parsed coverage data for {len(coverage_data)} contig name variants")
    return coverage_data


def parse_bin_files(bins_dir):
    """Parse bin files to get contig-to-bin mapping"""
    print(f"[INFO] Parsing bin files in: {bins_dir}")
    
    bin_contigs = {}
    
    # Find all .fa files
    fa_files = glob.glob(os.path.join(bins_dir, "*.fa"))
    
    if not fa_files:
        print(f"[ERROR] No .fa files found in {bins_dir}")
        return bin_contigs
    
    print(f"[INFO] Found {len(fa_files)} bin files")
    
    for fa_file in fa_files:
        bin_name = os.path.splitext(os.path.basename(fa_file))[0]
        contigs = []
        
        with open(fa_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Extract contig name (remove '>' and any additional info)
                    contig_name = line[1:].split()[0]
                    contigs.append(contig_name)
        
        bin_contigs[bin_name] = contigs
        print(f"[INFO] Bin {bin_name}: {len(contigs)} contigs")
    
    return bin_contigs


def calculate_bin_coverage(bin_contigs, coverage_data):
    """Calculate weighted average coverage for each bin"""
    print("[INFO] Calculating bin coverage...")
    
    bin_coverage = {}
    
    for bin_name, contigs in bin_contigs.items():
        total_weighted_coverage = 0
        total_length = 0
        
        for contig in contigs:
            if contig in coverage_data:
                contig_len = coverage_data[contig]['length']
                contig_cov = coverage_data[contig]['coverage']
                
                total_weighted_coverage += contig_cov * contig_len
                total_length += contig_len
            else:
                print(f"[WARNING] Contig {contig} from bin {bin_name} not found in coverage data")
        
        if total_length > 0:
            bin_coverage[bin_name] = total_weighted_coverage / total_length
        else:
            bin_coverage[bin_name] = 0.0
            print(f"[WARNING] No coverage data found for bin {bin_name}")
    
    print(f"[INFO] Calculated coverage for {len(bin_coverage)} bins")
    return bin_coverage


def run_checkm2(bins_dir, output_dir, threads, script_dir):
    """Run CheckM2 analysis using separate shell script"""
    print("[INFO] Running CheckM2 analysis...")
    
    checkm_output = os.path.join(output_dir, 'checkm2_output')
    # Don't remove here - let the shell script handle it
    
    # Path to the CheckM2 runner script
    checkm2_script = os.path.join(script_dir, 'run_checkm2.sh')
    
    # Make sure the script is executable
    run_command(f"chmod +x {checkm2_script}")
    
    # Run CheckM2 using the shell script
    cmd = f"bash {checkm2_script} {bins_dir} {checkm_output} {threads}"
    result = run_command(cmd, check=False)
    
    if result.returncode != 0:
        print("[ERROR] CheckM2 analysis failed")
        print(f"[ERROR] CheckM2 stderr: {result.stderr if result.stderr else 'No error message'}")
        return None
    
    # Parse CheckM2 results
    results_file = os.path.join(checkm_output, 'quality_report.tsv')
    if os.path.exists(results_file):
        df = pd.read_csv(results_file, sep='\t')
        print(f"[INFO] CheckM2 results loaded for {len(df)} bins")
        return df
    else:
        print("[ERROR] CheckM2 quality_report.tsv not found")
        return None


def parse_gtdbtk_results(gtdbtk_dir):
    """Parse GTDB-Tk classification results"""
    print(f"[INFO] Parsing GTDB-Tk results in: {gtdbtk_dir}")
    
    # Find GTDB-Tk summary files (both bacterial and archaeal)
    bac_files = glob.glob(os.path.join(gtdbtk_dir, "gtdbtk.bac120.summary.tsv"))
    ar_files = glob.glob(os.path.join(gtdbtk_dir, "gtdbtk.ar53.summary.tsv"))
    
    all_files = bac_files + ar_files
    
    if not all_files:
        print(f"[ERROR] No GTDB-Tk summary files found in {gtdbtk_dir}")
        return None
    
    print(f"[INFO] Found GTDB-Tk files: {[os.path.basename(f) for f in all_files]}")
    
    # Combine all results
    all_dfs = []
    
    for summary_file in all_files:
        print(f"[INFO] Reading GTDB-Tk results from: {summary_file}")
        try:
            df = pd.read_csv(summary_file, sep='\t')
            all_dfs.append(df)
            print(f"[INFO] Loaded {len(df)} bins from {os.path.basename(summary_file)}")
        except Exception as e:
            print(f"[ERROR] Failed to read GTDB-Tk file {summary_file}: {e}")
    
    if not all_dfs:
        return None
    
    # Combine all dataframes
    combined_df = pd.concat(all_dfs, ignore_index=True)
    print(f"[INFO] Total GTDB-Tk results: {len(combined_df)} bins")
    
    return combined_df


def split_classification(classification):
    """Split taxonomic classification into individual levels"""
    levels = {'d': '', 'p': '', 'c': '', 'o': '', 'f': '', 'g': '', 's': ''}
    
    if pd.isna(classification) or not classification:
        return levels
    
    # Split by semicolon
    parts = classification.split(';')
    
    for part in parts:
        part = part.strip()
        if '__' in part:
            level, name = part.split('__', 1)
            if level in levels:
                levels[level] = name
    
    return levels


def generate_annotation_metrics(gtdbtk_df, checkm_df, bin_coverage, output_dir):
    """Generate the final annotation metrics table"""
    print("[INFO] Generating annotation metrics table...")
    
    results = []
    
    # Get all unique bin names from all sources
    all_bins = set()
    
    if gtdbtk_df is not None:
        # Use names directly - no .fa extension to remove
        gtdbtk_df['bin_name'] = gtdbtk_df['user_genome']
        all_bins.update(gtdbtk_df['bin_name'].tolist())
    
    if checkm_df is not None:
        # Use names directly - no .fa extension to remove
        checkm_df['bin_name'] = checkm_df['Name']
        all_bins.update(checkm_df['bin_name'].tolist())
    
    all_bins.update(bin_coverage.keys())
    
    print(f"[INFO] Processing {len(all_bins)} unique bins")
    
    for bin_name in sorted(all_bins):
        result = {'bin_name': bin_name}
        
        # Add bin coverage
        result['bin_coverage'] = bin_coverage.get(bin_name, np.nan)
        
        # Add completeness and contamination from CheckM2
        if checkm_df is not None:
            checkm_row = checkm_df[checkm_df['bin_name'] == bin_name]
            if len(checkm_row) > 0:
                result['completeness'] = checkm_row.iloc[0]['Completeness']
                result['contamination'] = checkm_row.iloc[0]['Contamination']
            else:
                result['completeness'] = np.nan
                result['contamination'] = np.nan
        else:
            result['completeness'] = np.nan
            result['contamination'] = np.nan
        
        # Add taxonomic classification
        if gtdbtk_df is not None:
            gtdbtk_row = gtdbtk_df[gtdbtk_df['bin_name'] == bin_name]
            if len(gtdbtk_row) > 0:
                classification = gtdbtk_row.iloc[0]['classification']
                tax_levels = split_classification(classification)
                result.update(tax_levels)
            else:
                # Add empty taxonomic levels
                for level in ['d', 'p', 'c', 'o', 'f', 'g', 's']:
                    result[level] = ''
        else:
            # Add empty taxonomic levels
            for level in ['d', 'p', 'c', 'o', 'f', 'g', 's']:
                result[level] = ''
        
        results.append(result)
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Reorder columns
    column_order = ['bin_name', 'bin_coverage', 'completeness', 'contamination', 'd', 'p', 'c', 'o', 'f', 'g', 's']
    df = df[column_order]
    
    # Save results
    output_file = os.path.join(output_dir, 'annotation_metrics.tsv')
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"[INFO] Annotation metrics table saved to: {output_file}")
    print(f"[INFO] Total bins in table: {len(df)}")
    
    # Print summary
    print("\n[SUMMARY]")
    print(f"Bins with coverage data: {df['bin_coverage'].notna().sum()}")
    print(f"Bins with completeness data: {df['completeness'].notna().sum()}")
    print(f"Bins with contamination data: {df['contamination'].notna().sum()}")
    print(f"Bins with taxonomic classification: {(df['d'] != '').sum()}")
    
    return df


def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    print(f"[INFO] Starting annotation metrics generation...")
    print(f"[INFO] GTDB-Tk directory: {args.gtdbtk_dir}")
    print(f"[INFO] Bins directory: {args.bins_dir}")
    print(f"[INFO] Coverage file: {args.coverage_file}")
    print(f"[INFO] Output directory: {args.outdir}")
    print(f"[INFO] Threads: {args.threads}")
    
    # Check input files/directories
    if not os.path.exists(args.gtdbtk_dir):
        print(f"[ERROR] GTDB-Tk directory not found: {args.gtdbtk_dir}")
        sys.exit(1)
    
    if not os.path.exists(args.bins_dir):
        print(f"[ERROR] Bins directory not found: {args.bins_dir}")
        sys.exit(1)
    
    if not os.path.exists(args.coverage_file):
        print(f"[ERROR] Coverage file not found: {args.coverage_file}")
        sys.exit(1)
    
    # Parse coverage data
    coverage_data = parse_coverage_file(args.coverage_file)
    if not coverage_data:
        print("[ERROR] No coverage data could be parsed")
        sys.exit(1)
    
    # Parse bin files
    bin_contigs = parse_bin_files(args.bins_dir)
    if not bin_contigs:
        print("[ERROR] No bin files could be parsed")
        sys.exit(1)
    
    # Calculate bin coverage
    bin_coverage = calculate_bin_coverage(bin_contigs, coverage_data)
    
    # Run CheckM2
    print("\n[INFO] Running CheckM2...")
    checkm_df = run_checkm2(args.bins_dir, args.outdir, args.threads, args.script_dir)
    
    # Parse GTDB-Tk results
    print("\n[INFO] Parsing GTDB-Tk results...")
    gtdbtk_df = parse_gtdbtk_results(args.gtdbtk_dir)
    
    # Generate final metrics table
    print("\n[INFO] Generating final annotation metrics table...")
    df = generate_annotation_metrics(gtdbtk_df, checkm_df, bin_coverage, args.outdir)
    
    print("\n[INFO] Annotation metrics generation completed successfully!")
    print(f"[INFO] Results saved to: {args.outdir}/annotation_metrics.tsv")


if __name__ == '__main__':
    main()