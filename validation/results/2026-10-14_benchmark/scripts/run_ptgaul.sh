#!/usr/bin/env bash
# ptGAUL (bioconda, Flye 2.7 as pinned by its recipe) on the six Col-0 lines, R10.4.1, with the defaults but the
# plastome size (-g, NC_000932.1 length) and threads.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_ptgaul
B=${BENCH:?the benchmark folder}; W=$B/col0/compare/ptgaul; mkdir -p $W; cd $W
for f in $B/col0_lines/ont/*.fastq.gz; do s=$(basename $f | cut -d_ -f1,2)
  [ -e $s/done ] && continue; rm -rf $s; mkdir -p $s
  /usr/bin/time -f "%e s, %M kB" -o $s.time ptGAUL.sh -r $B/refs/NC_000932.1.fasta -l $f -t 32 -g 154478 -o $s > $s.log 2>&1 \
    && touch $s/done && echo "$s ok" || echo "$s FAILED"
done
echo DONE
