#!/usr/bin/env bash
# BACoN on the 32 human samples' chrM reads (1000 Genomes long-read consortium, R10.4.1; reads already limited to
# those aligned to chrM): templated and Flye, each compared with SKA2 (default) and minimap2. Reference: rCRS.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate BACoN_eval
B=${BENCH:?the benchmark folder}; W=$B/human_mt/bacon; cd $W; ref=$B/refs/NC_012920.1.gb
bacon --version
go() { run=$1; shift; echo "[$(date +%T)] $run $*"
  /usr/bin/time -a -f "%e s, %M kB: $*" -o $run.time bacon -r $ref -i reads -o $run "$@" -t 32 -p 8 >> $run.log 2>&1 && echo "  ok" || echo "  FAILED"; }
go samtools
go samtools --snp-method minimap2
go flye -a flye
go flye -a flye --snp-method minimap2
echo "[$(date +%T)] DONE"
