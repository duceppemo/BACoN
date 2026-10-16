#!/usr/bin/env bash
# Plastome truth of the 24 accessions: GetOrganelle 1.7.7.1 on the Illumina reads (embplant_pt, its recommended
# options for plant plastomes: -R 15 -k 21,45,65,85,105). Four accessions at a time, 12 threads each; starts after
# the Illumina mapping (illumina.sh).
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_getorganelle
B=${BENCH:?the benchmark folder}; W=$B/acc/truth/getorganelle; mkdir -p $W; cd $W
until grep -q DONE $B/acc/truth/illumina.log; do sleep 120; done
one() { r1=$1; s=$(basename $r1 | sed 's/_ERR[0-9]*_1.fastq.gz//'); r2=${r1/_1.fastq.gz/_2.fastq.gz}
  ls $s/*path_sequence.fasta > /dev/null 2>&1 && { echo "$s done"; return; }
  get_organelle_from_reads.py -1 $r1 -2 $r2 -o $s -F embplant_pt -R 15 -k 21,45,65,85,105 -t 12 --overwrite > $s.log 2>&1 \
    && echo "$s ok $(ls $s/*path_sequence.fasta | wc -l) path(s)" || echo "$s FAILED"; }
export -f one
ls $B/accessions/illumina/*_1.fastq.gz | xargs -P 4 -I{} bash -c 'cd '"$W"' && one {}'
echo DONE
