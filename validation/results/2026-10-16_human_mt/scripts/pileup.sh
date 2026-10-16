#!/usr/bin/env bash
# Base counts of every sample's Illumina chrM alignment (pileup_calls.py), 8 at a time.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_snp
H=${BENCH:?the benchmark folder}/human_mt; T=$H/truth
for b in $H/illumina/*.chrM.bam; do [ -f $b.bai ] || samtools index $b; done
ls $H/illumina/*.chrM.bam | xargs -P 8 -I{} sh -c 'python '"$T"'/pileup_calls.py {} '"$T"'/chrM.fasta '"$T"'/pileup/$(basename {} .chrM.bam).calls.tsv.gz'
echo DONE
