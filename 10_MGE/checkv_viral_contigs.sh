#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

INPUT_FASTA="$ROOT/hg/10_MGE/results/10_MGE/mge/mge_reports/virus_no_provirus.fna"
OUTPUT_DIR="$SCRIPT_DIR/checkv_output"
CHECKV="$ROOT/conda_envs/checkv_env/bin/checkv"
CHECKV_DB="$ROOT/databases/checkv_db/checkv-db-v1.5"
THREADS="${THREADS:-16}"

export PATH="$(dirname "$CHECKV"):$PATH"

"$CHECKV" end_to_end \
    "$INPUT_FASTA" \
    "$OUTPUT_DIR" \
    -t "$THREADS" \
    -d "$CHECKV_DB"
