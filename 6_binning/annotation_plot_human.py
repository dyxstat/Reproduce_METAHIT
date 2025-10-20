#!/usr/bin/env python3
import pandas as pd
import os
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for server
import matplotlib.pyplot as plt
import numpy as np

def load_gtdb_annotations(annotation_dirs):
    """Load GTDB-Tk annotation results from each binning tool"""
    annotations = {}
    
    for tool, ann_dir in annotation_dirs.items():
        annotations[tool] = {}
        
        # Load bacterial annotations
        bac_file = os.path.join(ann_dir, 'gtdbtk.bac120.summary.tsv')
        if os.path.exists(bac_file):
            bac_df = pd.read_csv(bac_file, sep='\t')
            for _, row in bac_df.iterrows():
                bin_id = row['user_genome']
                classification = row['classification']
                annotations[tool][bin_id] = classification
        
        # Load archaeal annotations
        ar_file = os.path.join(ann_dir, 'gtdbtk.ar53.summary.tsv')
        if os.path.exists(ar_file):
            ar_df = pd.read_csv(ar_file, sep='\t')
            for _, row in ar_df.iterrows():
                bin_id = row['user_genome']
                classification = row['classification']
                annotations[tool][bin_id] = classification
        
        print(f"{tool}: {len(annotations[tool])} annotated bins")
    
    return annotations

def load_high_quality_bins(stats_files, completeness_threshold=50, contamination_threshold=10):
    """Load high-quality bin IDs from stats files"""
    high_quality_bins = {}
    
    for tool, filepath in stats_files.items():
        if not os.path.exists(filepath):
            print(f"[WARNING] Missing stats file: {filepath}")
            high_quality_bins[tool] = set()
            continue
            
        df = pd.read_csv(filepath, sep="\t")
        if 'completeness' not in df.columns or 'contamination' not in df.columns:
            print(f"[ERROR] Missing required columns in: {filepath}")
            high_quality_bins[tool] = set()
            continue
            
        # Filter for high-quality bins
        high_quality = df[
            (df['completeness'] >= completeness_threshold) & 
            (df['contamination'] <= contamination_threshold)
        ]
        
        # Extract bin IDs
        bin_ids = set(high_quality.iloc[:, 0].astype(str))
        high_quality_bins[tool] = bin_ids
        print(f"{tool}: {len(bin_ids)} high-quality bins")
    
    return high_quality_bins

def parse_taxonomic_classification(classification):
    """Parse GTDB taxonomic classification string into levels"""
    if pd.isna(classification) or classification == '':
        return {'Species': None, 'Genus': None, 'Family': None, 'Order': None}
    
    # Initialize all levels as None
    taxonomy = {'Species': None, 'Genus': None, 'Family': None, 'Order': None}
    
    levels = classification.split(';')
    
    for level in levels:
        if level.startswith('s__'):
            taxonomy['Species'] = level[3:] if level != 's__' else None
        elif level.startswith('g__'):
            taxonomy['Genus'] = level[3:] if level != 'g__' else None
        elif level.startswith('f__'):
            taxonomy['Family'] = level[3:] if level != 'f__' else None
        elif level.startswith('o__'):
            taxonomy['Order'] = level[3:] if level != 'o__' else None
    
    return taxonomy

def get_taxonomic_sets(annotations, high_quality_bins):
    """Get taxonomic sets for each tool and taxonomic level"""
    tool_taxa = {}
    
    for tool in annotations.keys():
        tool_taxa[tool] = {
            'Species': set(),
            'Genus': set(), 
            'Family': set(),
            'Order': set()
        }
        
        # Only count high-quality bins
        hq_bins = high_quality_bins.get(tool, set())
        
        for bin_id, classification in annotations[tool].items():
            if bin_id in hq_bins:  # Only process high-quality bins
                taxonomy = parse_taxonomic_classification(classification)
                
                for level in ['Species', 'Genus', 'Family', 'Order']:
                    if taxonomy[level] is not None:
                        tool_taxa[tool][level].add(taxonomy[level])
    
    return tool_taxa

def calculate_pairwise_overlaps(tool_taxa, reference_tool, comparison_tool):
    """Calculate overlaps between reference tool and one comparison tool for each taxonomic level"""
    overlaps = {}
    
    for level in ['Species', 'Genus', 'Family', 'Order']:
        ref_taxa = tool_taxa[reference_tool][level]
        comp_taxa = tool_taxa[comparison_tool][level]
        
        # Calculate overlaps
        overlaps[level] = {
            f'{reference_tool}_only': len(ref_taxa - comp_taxa),
            f'{comparison_tool}_only': len(comp_taxa - ref_taxa), 
            'Both': len(ref_taxa & comp_taxa)
        }
    
    return overlaps

def create_comparison_plots(tool_taxa, output_file="annotation_plot_sheep.pdf"):
    """Create Figure 4b style comparison plots - 3 plots with split-top bars"""
    
    # Set up matplotlib fonts - use system default instead of Arial to avoid warnings
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': 30,          # Base font size for ticks/axis
        'axes.labelsize': 30,     # Axis labels
        'axes.titlesize': 35,     # Subplot titles  
        'xtick.labelsize': 30,    # X-axis tick labels
        'ytick.labelsize': 30,    # Y-axis tick labels
        'legend.fontsize': 30     # Legend font size
    })
    
    # Set up 3 subplots - MUCH LARGER figure for readability
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    
    # Define comparisons and colors - reordered: bin3C, MetaCC, ImputeCC
    reference_tool = 'METAHIT'
    comparisons = [
        ('bin3C', '#9b59b6'),     # Purple
        ('MetaCC', '#f39c12'),    # Orange
        ('ImputeCC', '#e74c3c')   # Red
    ]
    
    metahit_color = '#3498db'  # Blue (always for METAHIT)
    both_color = '#7f8c8d'     # Gray (always for Both)
    
    levels = ['Species', 'Genus', 'Family', 'Order']
    bar_width = 0.6
    all_overlaps = {}
    
    for i, (comp_tool, comp_color) in enumerate(comparisons):
        ax = axes[i]
        
        # Calculate overlaps for this comparison
        overlaps = calculate_pairwise_overlaps(tool_taxa, reference_tool, comp_tool)
        all_overlaps[comp_tool] = overlaps
        
        # Get data
        both_values = [overlaps[level]['Both'] for level in levels]
        metahit_only_values = [overlaps[level][f'{reference_tool}_only'] for level in levels]
        comp_only_values = [overlaps[level][f'{comp_tool}_only'] for level in levels]
        
        x_pos = np.arange(len(levels))
        split_width = bar_width / 2
        
        # Create bars with split-top design
        # Bottom: Gray (Both) - full width
        bars_both = ax.bar(x_pos, both_values, bar_width, 
                          color=both_color, label='Both', alpha=0.8, edgecolor='black')
        
        # Top left: Other tool only
        bars_comp = ax.bar(x_pos - split_width/2, comp_only_values, split_width, 
                          bottom=both_values, color=comp_color, label=f'{comp_tool} only', 
                          alpha=0.8, edgecolor='black')
        
        # Top right: METAHIT only  
        bars_metahit = ax.bar(x_pos + split_width/2, metahit_only_values, split_width,
                             bottom=both_values, color=metahit_color, label=f'{reference_tool} only', 
                             alpha=0.8, edgecolor='black')
        
        # Add value labels for bars
        for j, level in enumerate(levels):
            # Label for 'Both' values in center of gray bars
            if both_values[j] > 0:
                ax.text(x_pos[j], both_values[j]/2, str(both_values[j]),
                       ha='center', va='center', color='white', fontsize=25)
            
            # Label for comparison tool values above left bar
            if comp_only_values[j] > 0:
                y_pos = both_values[j] + comp_only_values[j] + 0.5
                ax.text(x_pos[j] - split_width/2, y_pos, str(comp_only_values[j]),
                       ha='center', va='bottom', color='black', fontsize=25)
            
            # Label for METAHIT values above right bar
            if metahit_only_values[j] > 0:
                y_pos = both_values[j] + metahit_only_values[j] + 0.5
                ax.text(x_pos[j] + split_width/2, y_pos, str(metahit_only_values[j]),
                       ha='center', va='bottom', color='black', fontsize=25)
        
        # Customize subplot - 25pt titles, no bold
        ax.set_ylabel('Number of taxa', fontsize=30)
        ax.set_title(f'{reference_tool} vs {comp_tool}', fontsize=35)  # No fontweight='bold'
        
        # Set x-axis
        ax.set_xticks(x_pos)
        ax.set_xticklabels(levels, fontsize=30)
        
        # Add grid
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_axisbelow(True)
        
        # Set y-axis
        max_height = max([both_values[j] + max(metahit_only_values[j], comp_only_values[j]) for j in range(len(levels))])
        ax.set_ylim(0, max_height * 1.1 if max_height > 0 else 10)
    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Annotation plot saved as '{output_file}'")
    
    return all_overlaps

def print_diversity_statistics(tool_taxa, all_overlaps):
    """Print detailed diversity statistics"""
    print("\n" + "="*60)
    print("TAXONOMIC DIVERSITY STATISTICS FOR sheep DATASET")
    print("="*60)
    
    # Individual tool diversity
    print("\nTaxonomic diversity per tool:")
    tools = ['METAHIT', 'ImputeCC', 'MetaCC', 'bin3C']
    levels = ['Species', 'Genus', 'Family', 'Order']
    
    for level in levels:
        print(f"\n{level} diversity:")
        for tool in tools:
            count = len(tool_taxa.get(tool, {}).get(level, set()))
            print(f"  {tool:12}: {count:3d} unique {level.lower()}")
    
    # Pairwise overlap analysis
    print(f"\nPairwise overlap analysis (METAHIT vs each tool):")
    for comp_tool in ['ImputeCC', 'bin3C', 'MetaCC']:
        print(f"\nMETAHIT vs {comp_tool}:")
        overlaps = all_overlaps[comp_tool]
        for level in levels:
            print(f"  {level}:")
            for category, count in overlaps[level].items():
                print(f"    {category:15}: {count:3d}")
    
    print("\n" + "="*60)

def main():
    # Configuration
    annotation_dirs = {
        "ImputeCC": "../annotation_binning/annotation_imputecc",
        "MetaCC": "../annotation_binning/annotation_metacc", 
        "bin3C": "../annotation_binning/annotation_bin3c",
        "METAHIT": "../annotation_binning/annotation_metahit"
    }
    
    stats_files = {
        "ImputeCC": "binsB.stats",
        "MetaCC": "binsA.stats",
        "bin3C": "binsC.stats", 
        "METAHIT": "metawrap_50_10_bins.stats"
    }
    
    print("Step 1: Loading GTDB-Tk annotations...")
    annotations = load_gtdb_annotations(annotation_dirs)
    
    print("\nStep 2: Loading high-quality bin information...")
    high_quality_bins = load_high_quality_bins(stats_files)
    
    print("\nStep 3: Extracting taxonomic sets...")
    tool_taxa = get_taxonomic_sets(annotations, high_quality_bins)
    
    print("\nStep 4: Creating annotation comparison plots...")
    current_dir = os.getcwd()
    output_path = os.path.join(current_dir, "annotation_plot_human.pdf")
    all_overlaps = create_comparison_plots(tool_taxa, output_path)
    
    print("\nStep 5: Printing statistics...")
    print_diversity_statistics(tool_taxa, all_overlaps)
    
    print(f"\nSUCCESS: Annotation plot saved as '{output_path}'")
    print(f"File location: {output_path}")
    print("Download this file to view the plot locally")

if __name__ == "__main__":
    main()