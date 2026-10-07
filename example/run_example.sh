#!/usr/bin/env bash
# Generate the simulated example (an annotated 30 kb circular reference, four samples with known SNPs and a sample
# metadata file), run BACoN on it with the default settings, and check the results against the truth: the SNP
# distances, the region, gene and effect of every SNP, and the report's annotation and metadata.
#
# Usage: bash example/run_example.sh [OUTPUT_FOLDER] [THREADS]
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
out=${1:-example_output}
threads=${2:-4}

python "$here/make_example.py" "$out/data"
bacon -r "$out/data/reference.fasta" -i "$out/data/reads" -o "$out/bacon" -t "$threads" -p 2 \
    --annotation "$out/data/reference.gb" --metadata "$out/data/metadata.tsv"

python "$here/check_example.py" "$out"
