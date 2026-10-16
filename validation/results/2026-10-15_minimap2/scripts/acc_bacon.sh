#!/usr/bin/env bash
# BACoN on the 24 accessions' R9.4.1 reads: plastome templated and Flye; mitogenome templated and Flye with
# --filter-by-reference (as recommended for mitogenomes and whole-genome reads).
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate BACoN_eval
B=${BENCH:?the benchmark folder}; W=$B/acc/bacon; cd $W
bacon --version
go() { run=$1; shift; echo "[$(date +%T)] $run"
  /usr/bin/time -f "%e s, %M kB" -o $run.time bacon -i reads -o $run "$@" -t 16 -p 2 > $run.log 2>&1 && echo "  ok" || echo "  FAILED"; }
go cp_samtools -r $B/refs/NC_000932.1.gb
go mt_samtools_fref -r $B/refs/NC_037304.1.gb --filter-by-reference
go cp_flye -r $B/refs/NC_000932.1.gb -a flye
go mt_flye_fref -r $B/refs/NC_037304.1.gb -a flye --filter-by-reference
echo "[$(date +%T)] DONE"
