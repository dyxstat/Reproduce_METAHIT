#!/usr/bin/env python3
"""Generate the statistics used as input to the bin-overlap Sankey plots."""

import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]

MASH = ROOT / "conda_envs/gtdbtk-2.4.0/bin/mash"
if not MASH.exists():
    found = shutil.which("mash")
    if not found:
        raise FileNotFoundError(
            "mash was not found. Expected repo-local path "
            f"{MASH} or an executable named mash in PATH."
        )
    MASH = Path(found)

COMPLETENESS_THRESHOLD = 50
CONTAMINATION_THRESHOLD = 10
MASH_DISTANCE_THRESHOLD = 0.01
SKETCH_SIZE = 10000

DATASETS = {
    "hg": {"label": "Human gut", "results": "results"},
    "sheep": {"label": "Sheep gut", "results": "results"},
    "pig": {"label": "Pig gut", "results": "results"},
    "cow": {"label": "Cow rumen", "results": "results"},
    "bovine": {"label": "Bovine skin", "results": "results"},
    "ww": {"label": "Wastewater", "results": "results"},
    "mat": {"label": "Hydrothermal mats", "results": "results"},
}

COMPARISONS = ["bin3C", "MetaCC", "ImputeCC"]


def dataset_paths(dataset, results_name):
    result_root = ROOT / dataset / "6_binning" / results_name
    metahict_root = result_root / "metahit"
    stats_files = {
        "METAHICT": metahict_root / "metahit_50_10_bins.stats",
        "bin3C": metahict_root / "work_files/binsC.stats",
        "MetaCC": metahict_root / "work_files/binsA.stats",
        "ImputeCC": metahict_root / "work_files/binsB.stats",
    }
    bin_dirs = {
        "METAHICT": metahict_root / "metahit_50_10_bins",
        "bin3C": result_root / "bin3c/fasta",
        "MetaCC": result_root / "metacc/BIN",
        "ImputeCC": result_root / "imputecc/FINAL_BIN",
    }
    return stats_files, bin_dirs


def load_medium_quality_bins(stats_files):
    bins = {}
    for tool, stats_file in stats_files.items():
        if not stats_file.exists():
            raise FileNotFoundError(f"Missing stats file for {tool}: {stats_file}")
        df = pd.read_csv(stats_file, sep="\t")
        required = {"completeness", "contamination"}
        if not required.issubset(df.columns):
            raise ValueError(f"Missing required columns in {stats_file}")
        medium_quality = df[
            (df["completeness"] >= COMPLETENESS_THRESHOLD)
            & (df["contamination"] < CONTAMINATION_THRESHOLD)
        ]
        bins[tool] = set(medium_quality.iloc[:, 0].astype(str))
    return bins


def get_bin_fasta_paths(medium_quality_bins, bin_dirs):
    paths = {}
    fasta_suffixes = {".fa", ".fasta", ".fna"}
    for tool, bin_ids in medium_quality_bins.items():
        bin_dir = bin_dirs[tool]
        if not bin_dir.exists():
            raise FileNotFoundError(f"Missing bin directory for {tool}: {bin_dir}")
        paths[tool] = {}
        for fasta in sorted(bin_dir.iterdir()):
            if fasta.suffix not in fasta_suffixes:
                continue
            bin_id = fasta.stem
            if bin_id in bin_ids:
                paths[tool][bin_id] = fasta
        missing = bin_ids - set(paths[tool])
        if missing:
            raise FileNotFoundError(
                f"{tool}: {len(missing)} medium-quality bins do not have FASTA files. "
                f"Examples: {', '.join(sorted(missing)[:5])}"
            )
    return paths


def sketch_bins(bin_paths, sketch_dir):
    sketch_paths = {}
    for tool, paths in bin_paths.items():
        prefix = sketch_dir / tool
        cmd = [
            str(MASH),
            "sketch",
            "-s",
            str(SKETCH_SIZE),
            "-o",
            str(prefix),
        ] + [str(path) for _bin_id, path in sorted(paths.items())]
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        sketch_paths[tool] = Path(f"{prefix}.msh")
    return sketch_paths


def bin_id_from_mash_name(name):
    return Path(name).stem


def mash_distances(sketch_a, sketch_b):
    cmd = [str(MASH), "dist", str(sketch_a), str(sketch_b)]
    result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        yield bin_id_from_mash_name(parts[0]), bin_id_from_mash_name(parts[1]), float(parts[2])


def calculate_pairwise_overlap(metahit_sketch, other_sketch):
    matched_metahit = set()
    matched_other = set()
    matched_pairs = []

    for metahit_bin, other_bin, distance in mash_distances(metahit_sketch, other_sketch):
        if distance < MASH_DISTANCE_THRESHOLD:
            matched_metahit.add(metahit_bin)
            matched_other.add(other_bin)
            matched_pairs.append((metahit_bin, other_bin, distance))

    return matched_metahit, matched_other, matched_pairs


def main():
    summary_lines = [
        "dataset\tcomparison\tmetahict_total\t"
        "binner_total\toverlap\tmetahict_only\tbinner_only"
    ]
    pair_lines = ["dataset\tcomparison\tmetahict_bin\tbinner_bin\tmash_distance"]

    for dataset, config in DATASETS.items():
        results_name = config["results"]
        stats_files, bin_dirs = dataset_paths(dataset, results_name)
        medium_quality_bins = load_medium_quality_bins(stats_files)
        bin_paths = get_bin_fasta_paths(medium_quality_bins, bin_dirs)

        with tempfile.TemporaryDirectory(prefix=f"sankey_{dataset}_mash_", dir="/tmp") as tmp:
            sketch_paths = sketch_bins(bin_paths, Path(tmp))
            for tool in COMPARISONS:
                print(f"Calculating {dataset}: METAHICT vs {tool} Mash overlaps...")
                matched_metahict, matched_other, matched_pairs = calculate_pairwise_overlap(
                    sketch_paths["METAHICT"], sketch_paths[tool]
                )
                metahict_total = len(medium_quality_bins["METAHICT"])
                binner_total = len(medium_quality_bins[tool])
                overlap = len(matched_other)
                metahict_only = metahict_total - len(matched_metahict)
                binner_only = binner_total - len(matched_other)
                summary_lines.append(
                    f"{dataset}\tMETAHICT_vs_{tool}\t"
                    f"{metahict_total}\t{binner_total}\t{overlap}\t{metahict_only}\t{binner_only}"
                )
                for metahict_bin, other_bin, distance in sorted(
                    matched_pairs, key=lambda item: (item[1], item[0], item[2])
                ):
                    pair_lines.append(
                        f"{dataset}\tMETAHICT_vs_{tool}\t{metahict_bin}\t{other_bin}\t{distance:.8f}"
                    )

    summary_file = SCRIPT_DIR / "sankey_stats.tsv"
    pairs_file = SCRIPT_DIR / "sankey_matched_pairs.tsv"
    summary_file.write_text("\n".join(summary_lines) + "\n")
    pairs_file.write_text("\n".join(pair_lines) + "\n")
    print(f"Saved: {summary_file}")
    print(f"Saved: {pairs_file}")


if __name__ == "__main__":
    main()
