#!/usr/bin/env bash
# TIPPo v2.4 on the six Col-0 lines, R10.4.1 reads (-p onthq), reference-free, organelle mode (plastome and
# mitogenome). Three lines at a time, 16 threads each.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_tipp
B=${BENCH:?the benchmark folder}; W=$B/col0/compare/tippo_ont; mkdir -p $W; cd $W
one() { f=$1; s=$(basename $f | cut -d_ -f1,2); mkdir -p $s; cd $s
  [ -e done ] || { /usr/bin/time -f "%e s, %M kB" -o ../$s.time TIPPo.v2.4.pl -f $f -p onthq -g organelle -t 16 > tippo.log 2>&1 && touch done; }
  echo "$s $( [ -e done ] && echo ok || echo FAILED)"; }
export -f one
ls $B/col0_lines/ont/*.fastq.gz | xargs -P 3 -I{} bash -c 'cd '"$W"' && one {}'
echo DONE
