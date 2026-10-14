#!/usr/bin/env bash
# BACoN on the six Col-0 MA32 lines, ONT R10.4.1 (PRJEB112076): plastome and mitogenome, templated and Flye, with
# the default filter and with --filter-by-reference. BENCH: the folder with refs/ and col0_lines/ont/; the runs go
# to the current folder. Each Flye run, and each --filter-by-reference run, starts from a copy of the templated
# run's folder: the baiting is reused (resume).
set -u
R=$BENCH/refs
mkdir -p reads; for f in $BENCH/col0_lines/ont/*.fastq.gz; do ln -sf $f reads/$(basename $f | cut -d_ -f1,2).fastq.gz; done
go() { run=$1; shift; echo "[$(date +%T)] $run"
  /usr/bin/time -f "%e s, %M kB" -o $run.time bacon -i reads -o $run "$@" -t 48 -p 6 > $run.log 2>&1 && echo "  ok" || echo "  FAILED $run"; }
for org in cp mt; do
  ref=$R/$([ $org = cp ] && echo NC_000932.1 || echo NC_037304.1).gb
  go ${org}_samtools -r $ref
  [ -d ${org}_flye ] || rsync -a --exclude 3_assembled --exclude 4_compared ${org}_samtools/ ${org}_flye/
  go ${org}_flye -r $ref -a flye
  for a in samtools flye; do
    [ -d ${org}_${a}_fref ] || rsync -a --exclude 2_filtered --exclude 3_assembled --exclude 4_compared ${org}_samtools/ ${org}_${a}_fref/
    go ${org}_${a}_fref -r $ref $([ $a = flye ] && echo -a flye) --filter-by-reference
  done
done
