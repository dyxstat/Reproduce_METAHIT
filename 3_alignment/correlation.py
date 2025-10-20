#!/usr/bin/env python3
import scipy.stats as stats

# For correlations involving chimeric fraction, exclude sheep gut
ratio_3d_filtered = [0.26, 0.29, 0.17, 0.47, 0.43]  # Remove sheep gut
info_fraction_filtered = [0.42, 0.32, 0.40, 0.56, 0.51]  # Remove sheep gut

print("CORRELATION ANALYSIS FOR ALL THREE METRICS")
print("="*60)

# 2. 3D Ratio vs Info Fraction (5 datasets, excluding sheep gut)
correlation_2, p_value_2 = stats.pearsonr(ratio_3d_filtered, info_fraction_filtered)
print(f"3D Ratio vs Info Fraction (n=5, excluding sheep gut):")
print(f"  Pearson correlation: {correlation_2:.4f}")
print(f"  P-value: {p_value_2:.4f}")
print()

print("="*60)
print("SUMMARY:")
print(f"3D vs Info:       r = {correlation_2:.3f}, p = {p_value_2:.3f}")