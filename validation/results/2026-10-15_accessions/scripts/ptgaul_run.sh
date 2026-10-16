#!/usr/bin/env bash
# ptGAUL 1.0.5 on the 24 accessions' R9.4.1 reads (plastome; the defaults but the plastome size and threads), three
# at a time with 16 threads; starts after the GetOrganelle truth.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_ptgaul
B=${BENCH:?the benchmark folder}; W=$B/acc/ptgaul; cd $W
until grep -q DONE $B/acc/truth/getorganelle.log; do sleep 120; done
one() { f=$1; s=$(basename $f | sed 's/_ERR[0-9]*.fastq.gz//'); [ -e $s/done ] && { echo "$s done"; return; }
  rm -rf $s; mkdir -p $s
  ptGAUL.sh -r $B/refs/NC_000932.1.fasta -l $f -t 16 -g 154478 -o $s > $s.log 2>&1 && touch $s/done && echo "$s ok" || echo "$s FAILED"; }
export -f one; export B
ls $B/accessions/ont/*.fastq.gz | xargs -P 3 -I{} bash -c 'cd '"$W"' && one {}'
echo DONE
