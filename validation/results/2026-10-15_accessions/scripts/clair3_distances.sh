#!/usr/bin/env bash
# Distances of the Clair3 pipeline: the consensus sequences (reference coordinates, N at low depth) and the
# reference as one alignment per organelle, snp-dists (differences between A, C, G, T only: positions with N in
# either sequence do not count). Starts after clair3/run.sh.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_snp
B=${BENCH:?the benchmark folder}; W=$B/acc/compare; cd $W
until grep -q DONE $B/acc/clair3/run.log; do sleep 120; done
for o in cp:NC_000932.1 mt:NC_037304.1; do n=${o%%:*}; r=${o#*:}
  { sed '1s/.*/>Reference/' $B/refs/$r.fasta; cat $B/acc/clair3/*/$n.fasta; } > clair3_$n.aln
  snp-dists -q clair3_$n.aln > distances/clair3_$n.tsv; done
echo DONE
