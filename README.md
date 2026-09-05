# Instruction to reproduce results in METAHICT paper

Scripts used to reproduce the statistical analyses and plots in the manuscript *METAHICT enables comprehensive and flexible genome-resolved microbiome analysis with metagenomic Hi-C* are organized below by METAHICT module.

## 6. Binning Module

- `binning_result.py` — Compares the numbers of MAGs recovered by METAHICT, bin3C, MetaCC, and ImputeCC across six completeness and contamination thresholds.
- `ensemble_methods.py` — Compares the numbers of MAGs recovered by METAHICT, Binning_refiner, and DAS Tool across six completeness and contamination thresholds.
- `non_hic_binners.py` — Compares the numbers of MAGs recovered by Hi-C-based and non-Hi-C binners across six completeness and contamination thresholds.
- `sankey_stats.py` — Calculates shared and method-specific medium-quality MAG counts between METAHICT and each Hi-C-based binner.

## 7. Reassembly Module

- `all_pies.py` — Compares viral and plasmid contigs recovered from residual assemblies and bin-specific reassemblies across the five short-read datasets.
- `bin_quality.py` — Compares MAG completeness and contamination before and after bin-specific reassembly across the five short-read datasets.
- `real_em_analysis.py` — Visualizes the fitted insert-size mixture model and read-selection cutoff for one short-read dataset.
- `reassembly_comparison.py` — Compares MAG completeness and contamination among original bins, shotgun-only reassemblies, and shotgun-plus-Hi-C reassemblies and performs paired statistical tests.
- `simulation_em_analysis.py` — Evaluates fitted insert-size mixture models and read-selection cutoffs for simulated shotgun and Hi-C reads.

## 9. Annotation Module

- `annotation_binning.py` — Compares the taxonomic breadth of medium-quality MAGs recovered by METAHICT and each Hi-C-based binner in the human gut, pig gut, and wastewater datasets.
- `supp_annotation_binning.py` — Compares the taxonomic breadth of medium-quality MAGs recovered by METAHICT and each Hi-C-based binner in the sheep gut, cow rumen, bovine skin, and hydrothermal mats datasets.
- `phylum_stacked.py` — Summarizes the phylum-level composition of medium-quality MAGs across seven datasets.
- `dissimilarity.py` — Calculates and plots pairwise Bray–Curtis dissimilarities in phylum-level MAG composition across seven datasets.

## 10. MGE Module

- `fae_mags_heatmap.py` — Generates a contact heatmap for all Faecalibacterium MAGs in the human gut dataset.
- `phylum_contact_bar.py` — Summarizes viral and plasmid contigs associated with host MAGs by host phylum in the human gut dataset.
- `checkv_viral_contigs.sh` — Runs CheckV to evaluate viral contigs identified by the MGE module in the human gut dataset.
