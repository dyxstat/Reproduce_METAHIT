#!/usr/bin/env python3
import os, sys, glob, csv, subprocess

def run_cmd(cmd):
    print(f"[CMD] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def fasta_count(fa):
    if not os.path.exists(fa):
        return 0
    with open(fa) as f:
        return sum(1 for line in f if line.startswith(">"))

def main(bin_dir, unmapped, outdir, db, threads):
    bins_dir = bin_dir

    # All genomad outputs stay inside outdir
    os.makedirs(outdir, exist_ok=True)
    result_dir = os.path.join(outdir, "genomad_results")
    os.makedirs(result_dir, exist_ok=True)

    combined = os.path.join(outdir, "bins_combined.fa")
    outfile  = os.path.join(outdir, "bins_unmapped_mge.csv")

    # 1. Combine bins into one fasta
    print(f"[INFO] Combining bins into {combined}")
    with open(combined, "w") as fout:
        for fa in glob.glob(os.path.join(bins_dir, "*.fa")):
            with open(fa) as fin:
                fout.write(fin.read())

    # 2. Run genomad on unmapped + combined bins
    run_cmd(
        f"genomad end-to-end --cleanup --splits 8 "
        f"--threads {threads} "
        f"{unmapped} {os.path.join(result_dir,'unmapped')} {db}"
    )
    run_cmd(
        f"genomad end-to-end --cleanup --splits 8 "
        f"--threads {threads} "
        f"{combined} {os.path.join(result_dir,'bins')} {db}"
    )

    # 3. Parse genomad results from *_summary FASTA files
    unmapped_summary = glob.glob(os.path.join(result_dir, "unmapped", "*_summary"))
    bins_summary     = glob.glob(os.path.join(result_dir, "bins", "*_summary"))

    if not unmapped_summary or not bins_summary:
        raise RuntimeError("Could not find genomad summary outputs")

    us = unmapped_summary[0]
    bs = bins_summary[0]

    unmapped_virus   = fasta_count(os.path.join(us, os.path.basename(us).replace("_summary","_virus.fna")))
    unmapped_plasmid = fasta_count(os.path.join(us, os.path.basename(us).replace("_summary","_plasmid.fna")))
    bins_virus       = fasta_count(os.path.join(bs, os.path.basename(bs).replace("_summary","_virus.fna")))
    bins_plasmid     = fasta_count(os.path.join(bs, os.path.basename(bs).replace("_summary","_plasmid.fna")))

    unmapped_mge = unmapped_virus + unmapped_plasmid
    bins_mge     = bins_virus + bins_plasmid

    # 4. Write output CSV
    with open(outfile, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["source", "num_mge_contigs"])
        writer.writerow(["unmapped_final_contigs", unmapped_mge])
        writer.writerow(["bins_combined", bins_mge])

    print(f"[INFO] unmapped: virus={unmapped_virus}, plasmid={unmapped_plasmid}, total={unmapped_mge}")
    print(f"[INFO] bins: virus={bins_virus}, plasmid={bins_plasmid}, total={bins_mge}")
    print(f"[INFO] Wrote summary to {outfile}")

if __name__ == "__main__":
    if len(sys.argv) < 11:
        print("Usage: python genomad_bins_unmapped.py "
              "--bin <BIN> --unmapped <UNMAPPED> --outdir <OUTDIR> --db <DB> -t <THREADS>")
        sys.exit(1)

    args = sys.argv
    bin_dir  = args[args.index("--bin")+1]
    unmapped = args[args.index("--unmapped")+1]
    outdir   = args[args.index("--outdir")+1]
    db       = args[args.index("--db")+1]
    threads  = int(args[args.index("-t")+1])

    main(bin_dir, unmapped, outdir, db, threads)
