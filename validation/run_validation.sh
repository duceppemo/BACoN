#!/usr/bin/env bash
# Run BACoN on the simulated scenarios with every assembler and SNP method, then score against the truth.
#
# Usage: bash validation/run_validation.sh OUTPUT_FOLDER [THREADS]
#
# Run it in the BACoN conda environment (environment.yml). The simulated reads are generated first (fixed seed:
# the same files every time); each run writes into OUTPUT_FOLDER/<scenario>/<run> and the scores go to
# OUTPUT_FOLDER/SCORES_<scenario>.md and .json. See validation/README.md.
set -euo pipefail

out=$(realpath -m "${1:?output folder}")
threads=${2:-16}
here=$(cd "$(dirname "$0")" && pwd)

if [ -d "$out/plastid" ]; then  # Resumed runs would skip steps and make the timings meaningless
    echo "error: $out already has results; use a new folder (simulated reads can be reused: copy $out/sim)" >&2
    exit 1
fi
mkdir -p "$out"
if [ ! -f "$out/sim/rdna/truth.json" ]; then
    python "$here/simulate.py" "$out/sim"
fi

bacon_run() {  # bacon_run SCENARIO RUN [bacon options...]
    local scenario=$1 run=$2
    shift 2
    local dir="$out/$scenario/$run"
    mkdir -p "$out/$scenario"
    echo "[$(date +%T)] $scenario / $run: bacon $*"
    # One line per bacon call, appended: the first is the full run, later ones only the steps that changed.
    /usr/bin/time -a -f "%e s, %M kB: bacon $*" -o "$dir.time" \
        bacon -r "$out/sim/$scenario/reference.fasta" -i "$out/sim/$scenario/reads" -o "$dir" \
        -t "$threads" -p 4 "$@" >> "$dir.log" 2>&1 || echo "    FAILED (see $dir.log)"
}

parsnp_default() {  # Parsnp without -c (its default genome filtering), on existing assemblies
    local scenario=$1 dir="$out/$1/flye/4_compared/parsnp_default"
    echo "[$(date +%T)] $scenario / flye: parsnp without -c"
    rm -rf "$dir"
    mkdir -p "$(dirname "$dir")"
    local asm=("$out/$scenario/flye/3_assembled/all_assemblies/"*.fasta)
    parsnp -r "$out/$scenario/flye/reference.fasta" -d "${asm[@]}" -o "$dir" -p "$threads" > "$dir.log" 2>&1 || true
    if [ -f "$dir/parsnp.ggr" ]; then
        harvesttools -i "$dir/parsnp.ggr" -S "$dir/snps.raw.fasta" >> "$dir.log" 2>&1
        python - "$dir" "$out/$scenario/flye" <<'PY'
import sys
from pathlib import Path
from bacon.compare import clean_alignment, snp_distances, write_distances
d, run = Path(sys.argv[1]), Path(sys.argv[2])
rename = {p.name: p.stem for p in (run / "3_assembled" / "all_assemblies").glob("*.fasta")}
rename.update({"reference.fasta.ref": "Reference", "reference.fasta": "Reference"})
recs = clean_alignment(d / "snps.raw.fasta", d / "parsnp.snps.fasta", rename)
write_distances(d / "snp_distances.tsv", *snp_distances(recs))
PY
    fi
}

{
    echo "# Validation run $(date -Iseconds)"
    echo "bacon: $(bacon --version)"
    for tool in minimap2 filtlong flye myloasm samtools ska harvesttools; do
        echo "$tool: $("$tool" --version 2>&1 | head -1)"
    done
    echo "parsnp: $(parsnp --version 2>&1 | grep -o 'Parsnp [0-9.]*' | head -1)"
    echo "FastTree: $(FastTree -expert 2>&1 | grep -o 'FastTree [0-9.]*' | head -1)"
} > "$out/versions.txt"

for scenario in plastid linear rdna; do
    bacon_run "$scenario" samtools                  # the default: templated assembly, SKA2 core SNPs
    bacon_run "$scenario" flye -a flye
    bacon_run "$scenario" myloasm -a myloasm
    # Other comparisons of the Flye assemblies (the assembly step is reused from the checkpoint).
    bacon_run "$scenario" flye -a flye --snp-method parsnp
    bacon_run "$scenario" flye -a flye --ska-min-freq 0.5
    parsnp_default "$scenario"
    python "$here/score.py" "$out/sim/$scenario" "$out/$scenario"/*/ --markdown "$out/SCORES_$scenario.md" \
        --json "$out/SCORES_$scenario.json" > /dev/null
    echo "Scores: $out/SCORES_$scenario.md"
done
