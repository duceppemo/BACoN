#!/usr/bin/env bash
# Generate the simulated example (a 30 kb circular reference and four samples with known SNPs), run BACoN on it
# with the default settings, and check the SNP distances against the truth.
#
# Usage: bash example/run_example.sh [OUTPUT_FOLDER] [THREADS]
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
out=${1:-example_output}
threads=${2:-4}

python "$here/make_example.py" "$out/data"
bacon -r "$out/data/reference.fasta" -i "$out/data/reads" -o "$out/bacon" -t "$threads" -p 2

python - "$out/bacon" "$here" <<'PY'
import sys
from pathlib import Path

sys.path.insert(0, sys.argv[2])
from make_example import EXPECTED_DISTANCES

found = sorted(Path(sys.argv[1], "4_compared").glob("*/snp_distances.tsv"))
if not found:
    sys.exit("FAILED: no snp_distances.tsv")
lines = found[0].read_text().splitlines()
names = lines[0].split("\t")[1:]
dist = {row.split("\t")[0]: dict(zip(names, map(int, row.split("\t")[1:]))) for row in lines[1:]}
wrong = [f"{a}-{b}: {dist[a][b]} (expected {d})" for (a, b), d in EXPECTED_DISTANCES.items() if dist[a][b] != d]
if wrong:
    sys.exit("FAILED: " + "; ".join(wrong))
print(f"OK: all {len(EXPECTED_DISTANCES)} pairwise SNP distances match the truth ({found[0]})")
PY
