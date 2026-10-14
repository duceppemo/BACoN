#!/usr/bin/env bash
# Fresh BACoN runs for timing (the validation runs were resumed): six Col-0 lines, one at a time (-p 1), 16
# threads (as TIPPo, Oatk and CLAW per job). Plastome templated and Flye; mitogenome templated and Flye with
# --filter-by-reference. Waits for the Oatk runs to finish (an idle machine).
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate BACoN_eval
B=${BENCH:?the benchmark folder}; W=$B/col0/compare/bacon_timing; mkdir -p $W; cd $W
until grep -q DONE $B/col0/compare/oatk/run.log; do sleep 120; done
mkdir -p reads; for f in $B/col0_lines/ont/*.fastq.gz; do ln -sf $f reads/$(basename $f | cut -d_ -f1,2).fastq.gz; done
bacon --version
go() { run=$1; shift; echo "[$(date +%T)] $run"
  /usr/bin/time -f "%e s, %M kB" -o $run.time bacon -i reads -o $run "$@" -t 16 -p 1 > $run.log 2>&1 && echo "  ok" || echo "  FAILED"; }
go cp_samtools -r $B/refs/NC_000932.1.gb
go cp_flye -r $B/refs/NC_000932.1.gb -a flye
go mt_samtools_fref -r $B/refs/NC_037304.1.gb --filter-by-reference
go mt_flye_fref -r $B/refs/NC_037304.1.gb -a flye --filter-by-reference
echo "[$(date +%T)] DONE"
