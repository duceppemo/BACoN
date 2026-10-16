#!/usr/bin/env bash
# Base counts of every accession's Illumina alignment (acc_truth_illumina.sh) at each organelle position
# (acc_pileup_calls.py), then the truth tables (acc_truth.py). The reference: TAIR10.1 with the interval of
# nuc_mask.bed (the nuclear copy of the mitogenome on chromosome 2) replaced by N.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_snp
B=${BENCH:?the benchmark folder}; T=$B/acc/truth; S=$(cd "$(dirname "$0")" && pwd)
mkdir -p $T/pileup
ls $T/illumina/*.bam | xargs -P 8 -I{} sh -c \
    'python '"$S"'/acc_pileup_calls.py {} '"$B"'/refs/TAIR10.1.orgmask.fasta '"$T"'/pileup/$(basename {} .bam).calls.tsv.gz'
BENCH=$B python $S/acc_truth.py
