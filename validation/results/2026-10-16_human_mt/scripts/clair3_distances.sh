#!/usr/bin/env bash
# Distances of the Clair3 baseline: the consensus sequences (rCRS coordinates, N at low depth) and rCRS as one
# alignment, snp-dists (differences between A, C, G, T only). After clair3/run.sh.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_snp
H=${BENCH:?the benchmark folder}/human_mt; cd $H/compare
{ sed '1s/.*/>Reference/' $H/truth/chrM.fasta; cat $H/clair3/*/chrM.fasta; } > clair3.aln
snp-dists -q clair3.aln > clair3.distances.tsv
