#!/usr/bin/env bash
# Score every run (BACoN's validation runs of the 2026-10-13_mito record and compare/collect/*) against the HiFi
# truths with score.py; then the SKA2 distances (distances.py, in BACoN's environment).
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_truth
export BENCH=${BENCH:?the benchmark folder}; S=$(cd "$(dirname "$0")" && pwd); C=$BENCH/col0/compare; cd $C; mkdir -p scores
runs=$(ls -d $BENCH/col0/bacon/{cp,mt}_{samtools,flye}{,_fref} $C/collect/*)
for r in $runs; do
  n=$(basename $r); case $n in cp_*|mt_*) n=bacon_$n;; esac; org=$(case $n in *_cp|*cp_*) echo cp;; *) echo mt;; esac)
  python $S/score.py $r $org > scores/$n.tsv && echo "$n scored"
done
conda activate BACoN_eval; python $S/distances.py $runs
