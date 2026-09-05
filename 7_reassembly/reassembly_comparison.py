import re
import warnings
from decimal import Decimal
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import rankdata, wilcoxon


FILES = {
    "Original": Path("pre_reassembly.tsv"),
    "SG-only": Path("sg_only.tsv"),
    "SG+Hi-C": Path("metahict_reassembled.tsv"),
}
DATASET_NAME = "Dataset"
OUTPUT_DIR = Path("reassembly_results")
FIGURE_DPI = 300

METHODS = ["Original", "SG-only", "SG+Hi-C"]
PREFIXES = {"Original": "Original", "SG-only": "SG_only", "SG+Hi-C": "SG_HiC"}
COMPARISONS = [
    ("Original", "SG-only", "Reassembly-only effect"),
    ("SG-only", "SG+Hi-C", "Recovered-Hi-C incremental effect"),
    ("Original", "SG+Hi-C", "Total reassembly effect"),
]


def read_quality_file(path):
    separator = "," if ".csv" in Path(path).suffixes else "\t"
    frame = pd.read_csv(path, sep=separator, compression="infer", dtype=str, skip_blank_lines=False)
    columns = {re.sub(r"[^a-z0-9]", "", str(c).lower()): c for c in frame.columns}
    needed = [columns.get("completeness") or columns.get("completion"), columns.get("contamination")]
    if any(c is None for c in needed):
        raise ValueError(f"{path}: required columns are Completeness and Contamination.")
    out = pd.DataFrame({
        metric: pd.to_numeric(frame[column], errors="coerce")
        for metric, column in zip(("Completeness", "Contamination"), needed)
    })
    if not np.isfinite(out.to_numpy(dtype=float)).all():
        raise ValueError(f"{path}: completeness/contamination contains missing, nonnumeric or infinite values.")
    dataset_column = columns.get("dataset")
    out["Dataset"] = frame[dataset_column] if dataset_column else DATASET_NAME
    if out["Dataset"].isna().any() or out["Dataset"].astype(str).str.strip().eq("").any():
        raise ValueError(f"{path}: dataset labels cannot be missing.")
    out["Dataset"] = out["Dataset"].astype(str).str.strip()
    out["Bin"] = np.arange(1, len(out) + 1)
    return out


def bh_adjust(p_values):
    values = np.asarray(p_values, dtype=float)
    result = np.full(values.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(values))
    if not len(valid):
        return result
    if np.any((values[valid] < 0) | (values[valid] > 1)):
        raise ValueError("P-values must be between zero and one.")
    order = valid[np.argsort(values[valid], kind="stable")]
    adjusted = values[order] * len(values) / np.arange(1, len(order) + 1)
    result[order] = np.clip(np.minimum.accumulate(adjusted[::-1])[::-1], 0, 1)
    return result


def compute_statistics(binwise, datasets, METHODS, COMPARISONS):
    prefixes = {"Original": "Original", "SG-only": "SG_only", "SG+Hi-C": "SG_HiC"}
    metrics = ("Completeness", "Contamination")
    datasets, comparisons = list(datasets), list(COMPARISONS)
    if not datasets or not comparisons:
        raise ValueError("Provide at least one dataset and one comparison.")
    required = [f"{prefixes[m]}_{metric}" for m in METHODS for metric in metrics]
    missing = set(["Dataset", *required]) - set(binwise.columns)
    if missing:
        raise ValueError(f"Missing paired-data columns: {sorted(missing)}")
    if not np.isfinite(binwise[required].to_numpy(dtype=float)).all():
        raise ValueError("Statistical input must be the finite, complete paired cohort.")
    rows = []
    for dataset in datasets:
        subset = binwise.loc[binwise["Dataset"] == dataset]
        for metric in metrics:
            for before_method, after_method, comparison in comparisons:
                before = subset[f"{prefixes[before_method]}_{metric}"].to_numpy(dtype=float)
                after = subset[f"{prefixes[after_method]}_{metric}"].to_numpy(dtype=float)

                delta = np.array([
                    float(Decimal(str(a)) - Decimal(str(b)))
                    for a, b in zip(after, before)
                ])
                if not np.isfinite(delta).all():
                    raise ValueError("Paired subtraction produced nonfinite differences.")
                n = len(delta)
                improved = delta < 0 if metric == "Contamination" else delta > 0
                worsened = delta > 0 if metric == "Contamination" else delta < 0
                tied = delta == 0
                ranks = rankdata(np.abs(delta), method="average")
                positive, negative = ranks[delta > 0].sum(), ranks[delta < 0].sum()
                effect = ((positive - negative) / (positive + negative)
                          if positive + negative else (0.0 if n else np.nan))
                p_value, statistic, note = np.nan, np.nan, ""
                if n < 2:
                    status = "insufficient_pairs"
                    note = "At least two matched pairs are required."
                elif tied.all():
                    p_value, statistic, status = 1.0, 0.0, "all_zero_differences"
                else:
                    try:
                        with warnings.catch_warnings(record=True) as caught:
                            warnings.simplefilter("always")
                            test = wilcoxon(delta, alternative="two-sided", zero_method="pratt", method="auto")
                        p_value, statistic = float(test.pvalue), float(test.statistic)
                        note = "; ".join(str(item.message) for item in caught)
                        status = "ok_with_warning" if caught else "ok"
                        if not np.isfinite(p_value) or not 0 <= p_value <= 1:
                            p_value, status = np.nan, "invalid_test_result"
                            note = (note + "; " if note else "") + "Wilcoxon returned an invalid p-value."
                    except (ValueError, FloatingPointError, OverflowError) as exc:
                        status, note = "test_error", str(exc)
                rows.append({
                    "Dataset": dataset, "Metric": metric, "Comparison": comparison,
                    "BeforeMethod": before_method, "AfterMethod": after_method, "N": n,
                    "MedianBefore": float(np.median(before)) if n else np.nan,
                    "MedianAfter": float(np.median(after)) if n else np.nan,
                    "MedianDeltaAfterMinusBefore": float(np.median(delta)) if n else np.nan,
                    "MeanDeltaAfterMinusBefore": float(np.mean(delta)) if n else np.nan,
                    "ImprovedN": int(improved.sum()), "WorsenedN": int(worsened.sum()),
                    "TiedN": int(tied.sum()), "NonzeroN": int((~tied).sum()),
                    "ImprovedFraction": float(improved.mean()) if n else np.nan,
                    "WorsenedFraction": float(worsened.mean()) if n else np.nan,
                    "WilcoxonP": p_value, "WilcoxonStatistic": statistic,
                    "PairedRankBiserial": float(effect),
                    "EffectRankConvention": "Pratt; positive means after minus before > 0",
                    "Test": "Wilcoxon two-sided; zero_method=pratt; method=auto",
                    "TestStatus": status, "TestNote": note,
                })
    statistics = pd.DataFrame(rows)
    statistics["BH_Q"] = np.nan
    for _, family in statistics.groupby("Comparison", sort=False):
        p_values = family["WilcoxonP"].to_numpy(dtype=float)
        statistics.loc[family.index, "BH_Q"] = bh_adjust(p_values)
        statistics.loc[family.index, "BH_FamilyN"] = len(family)
        statistics.loc[family.index, "BH_TestedN"] = int(np.isfinite(p_values).sum())
    statistics[["BH_FamilyN", "BH_TestedN"]] = statistics[["BH_FamilyN", "BH_TestedN"]].astype(int)
    statistics["BH_Family"] = "Within comparison: Completeness and Contamination across all datasets"
    return statistics


def plot_reassembly_boxplot(long_df, stats_df, datasets, output_prefix, dpi=300):
    methods = ("Original", "SG-only", "SG+Hi-C")
    labels = ("Pre-reassembly", "SG-only", "METAHICT Reassembled")
    comparisons = (
        ("Original", "SG-only", "Reassembly-only effect"),
        ("SG-only", "SG+Hi-C", "Recovered-Hi-C incremental effect"),
        ("Original", "SG+Hi-C", "Total reassembly effect"),
    )
    palettes = {
        "Contamination": ("#cce6ff", "#66b3ff", "#0066cc"),
        "Completeness": ("#ffcccc", "#ff6666", "#cc0000"),
    }
    datasets = list(datasets)
    if not datasets or len(set(datasets)) != len(datasets):
        raise ValueError("datasets must contain at least one unique dataset name.")
    required = {"Dataset", "Bin", "Method", "Completeness", "Contamination"}
    if not required.issubset(long_df.columns):
        raise ValueError(f"Plot data are missing columns: {sorted(required - set(long_df.columns))}")
    required_stats = {"Dataset", "Metric", "Comparison", "BeforeMethod", "AfterMethod", "BH_Q"}
    if not required_stats.issubset(stats_df.columns):
        raise ValueError(f"Statistics are missing columns: {sorted(required_stats - set(stats_df.columns))}")
    if isinstance(dpi, bool) or not isinstance(dpi, (int, np.integer)) or dpi <= 0:
        raise ValueError("dpi must be a positive integer.")
    arrays = {metric: {method: [] for method in methods} for metric in palettes}
    for dataset in datasets:
        subset = long_df.loc[long_df["Dataset"] == dataset]
        if subset[["Bin", "Method"]].isna().any().any() or subset.duplicated(["Bin", "Method"]).any():
            raise ValueError(f"Dataset {dataset!r} has missing or duplicate bin/method identifiers.")
        cohorts = []
        for method in methods:
            sub = subset.loc[subset["Method"] == method]
            cohorts.append(set(sub["Bin"]))
            for metric in palettes:
                values = sub[metric].to_numpy(dtype=float)
                if not len(values) or not np.isfinite(values).all():
                    raise ValueError(f"Dataset {dataset!r}, {method}, {metric}: expected finite matched observations.")
                arrays[metric][method].append(values)
        if not all(cohort == cohorts[0] for cohort in cohorts[1:]):
            raise ValueError(f"Dataset {dataset!r} does not use the same bins for all three methods.")

    prefix = Path(output_prefix).expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    paths = {suffix: str(prefix) + f".{suffix}" for suffix in ("pdf", "png")}
    centers = np.arange(len(datasets), dtype=float) * 2.4
    positions = {method: centers + shift for method, shift in zip(methods, (-0.42, 0.0, 0.42))}
    width = max(11.0, 4.0 + 2.8 * len(datasets))
    rc = {"font.family": "DejaVu Sans", "font.size": 13, "axes.titlesize": 20,
          "axes.labelsize": 15, "xtick.labelsize": 13, "ytick.labelsize": 13,
          "pdf.fonttype": 42, "ps.fonttype": 42}
    with plt.rc_context(rc):
        fig, axes = plt.subplots(1, 2, figsize=(width, 7.2))
        try:
            for ax, (metric, colors) in zip(axes, palettes.items()):
                boxplots = {}
                for method, color in zip(methods, colors):
                    bp = ax.boxplot(arrays[metric][method], positions=positions[method], widths=0.32,
                                    patch_artist=True, showfliers=True, whis=1.5, manage_ticks=False,
                                    flierprops={"marker": ".", "markersize": 3, "alpha": 0.65})
                    for box in bp["boxes"]:
                        box.set(facecolor=color, edgecolor="#333333", linewidth=1.0)
                    for line in bp["whiskers"] + bp["caps"]:
                        line.set(color="#333333", linewidth=1.0)
                    for median in bp["medians"]:
                        median.set(color="#f4a261", linewidth=1.3)
                    boxplots[method] = bp
                visible = np.concatenate([values for method in methods for values in arrays[metric][method]])
                low, high = float(visible.min()), float(visible.max())
                scale = max(high - low, abs(high) * 0.05, 1.0)
                lower = min(0.0, low - 0.04 * scale) if metric == "Contamination" else low - 0.06 * scale

                for i, dataset in enumerate(datasets):
                    local_top = max(float(np.max(arrays[metric][m][i])) for m in methods)
                    for level, (before, after, comparison) in enumerate(comparisons):
                        rows = stats_df.loc[
                            (stats_df["Dataset"] == dataset) & (stats_df["Metric"] == metric)
                            & (stats_df["Comparison"] == comparison)
                            & (stats_df["BeforeMethod"] == before) & (stats_df["AfterMethod"] == after)
                        ]
                        if len(rows) != 1:
                            raise ValueError(f"Expected one statistics row for {dataset!r}, {metric}, {comparison}.")
                        q = float(rows.iloc[0]["BH_Q"])
                        if np.isfinite(q) and not 0 <= q <= 1:
                            raise ValueError("BH_Q values must be between zero and one, or missing.")
                        label = "NA" if not np.isfinite(q) else ("*" if q <= 0.05 else "NS")
                        base = local_top + (0.04 + 0.095 * level) * scale
                        top = base + 0.022 * scale
                        x1, x2 = positions[before][i], positions[after][i]
                        ax.plot([x1, x1, x2, x2], [base, top, top, base], color="black", linewidth=1.05)
                        ax.text((x1 + x2) / 2, top + 0.009 * scale, label, ha="center", va="bottom", fontsize=13)
                ax.set(title=metric, ylabel=metric, xlim=(centers[0] - 1.0, centers[-1] + 1.0),
                       ylim=(lower, high + 0.36 * scale))
                ax.set_xticks(centers, datasets)
                ax.legend(handles=[Patch(facecolor=color, edgecolor="#333333", label=label)
                                   for color, label in zip(colors, labels)],
                          loc="upper center", bbox_to_anchor=(0.5, -0.13), frameon=False, fontsize=9, ncol=3, handlelength=1.8, columnspacing=1.0)
            fig.subplots_adjust(left=0.075, right=0.985, top=0.90, bottom=0.20, wspace=0.25)
            for suffix, path in paths.items():
                fig.savefig(path, dpi=dpi, bbox_inches="tight")
        finally:
            plt.close(fig)
    return paths


frames = {method: read_quality_file(path) for method, path in FILES.items()}
reference = frames["Original"]["Dataset"].tolist()
if any(frame["Dataset"].tolist() != reference for frame in frames.values()):
    raise ValueError("All three files must have equal row counts and matching dataset order.")
if not reference:
    raise ValueError("The input files contain no data rows.")


wide = frames["Original"][["Dataset", "Bin"]].copy()
for method in METHODS:
    for metric in ("Completeness", "Contamination"):
        wide[f"{PREFIXES[method]}_{metric}"] = frames[method][metric].to_numpy(copy=True)

dataset_order = ["Human", "Pig", "Bovine", "Wastewater", "Mats"]
present = list(wide["Dataset"].drop_duplicates())
datasets = [d for d in dataset_order if d in present] + [d for d in present if d not in dataset_order]
for dataset in datasets:
    print(f"{dataset}: all {int(wide['Dataset'].eq(dataset).sum())} rows per method included.")

long_frames = []
for method in METHODS:
    prefix = PREFIXES[method]
    temp = wide[["Dataset", "Bin", f"{prefix}_Completeness", f"{prefix}_Contamination"]].rename(columns={f"{prefix}_Completeness": "Completeness", f"{prefix}_Contamination": "Contamination"})
    temp["Method"] = method
    long_frames.append(temp)
long_df = pd.concat(long_frames, ignore_index=True)

stats_df = compute_statistics(wide, datasets, METHODS, COMPARISONS)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
stats_df.to_csv(OUTPUT_DIR / "reassembly_statistics.tsv", sep="\t", index=False)
wide.to_csv(OUTPUT_DIR / "reassembly_matched_bins.tsv", sep="\t", index=False)
plot_reassembly_boxplot(long_df, stats_df, datasets, OUTPUT_DIR / "reassembly_comparison", dpi=FIGURE_DPI)
print(f"Saved figure and paired statistics in {OUTPUT_DIR.resolve()}")
