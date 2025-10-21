# Instruction to reproduce results in METAHIT paper
Scripts to generate the intermediate data and plots to illustrate the results of some METAHIT modules are available in their corresponding folders.

## 6. Binning Module  
The folder `6_binning` contains three scripts for evaluating the results of the binning module:

- **annotation_plot_human.py** — Loads GTDB-Tk annotations and high-quality bin statistics, compares taxonomic diversity across tools, and generates annotation overlap plots.  
- **count_methods_bins.py** — Counts the number of bins that meet completeness and contamination thresholds for each binning method.  
- **binning_result.py** — Plots the number of bins as a stacked bar chart.  

## 7. Reassembly Module  
The folder `7_reassembly` contains six scripts for evaluating the results of the reassembly module:

- **avg_comp.py** — Calculates the average completeness of high-quality bins before and after reassembly.  
- **avg_cont.py** — Calculates the average contamination of high-quality bins before and after reassembly.  
- **bin_quality.py** — Plots completeness and contamination curves before and after reassembly.  
- **cont_drop.py** — Compares bin-by-bin contamination levels between original and reassembled bins.  
- **genomad_bins_unmapped.py** — Runs geNomad on both unmapped contigs and reassembled bins to identify viral and plasmid contigs and outputs a summary table.  
- **reassembly_metrics.py** — Evaluates reassembly effectiveness across seven population-level metrics.

## 9. Annotation Module  
The folder `9_annotation` contains one script for evaluating the results of the annotation module:

- **annotation_metrics.py** — Combines per-bin coverage, completeness, contamination, and GTDB-Tk taxonomic classification into a table.

## 10. MGE Module  
The folder `10_MGE` contains three scripts for evaluating the results of the MGE module:

- **run_virgo_subset.py** — Runs Virgo per-contig on QC-passed viral contigs linked to selected bins and annotates them by taxonomic order.  
- **phylum_bin_contact.py** — Counts viral and plasmid contigs contacting bins across the top four phyla.  
- **phylum_contact_bar.py** — Plots the viral and plasmid contact counts as a stacked bar chart.

