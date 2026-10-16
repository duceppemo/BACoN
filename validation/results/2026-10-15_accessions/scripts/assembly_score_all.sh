#!/usr/bin/env bash
# Plastome assemblies of every pipeline against the GetOrganelle truths (assembly_score.py): BACoN templated and
# Flye, ptGAUL (oriented), the Clair3 consensus. Writes compare/scores_assembly_<set>.tsv.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_truth
export BENCH=${BENCH:?the benchmark folder}; A=$BENCH/acc; cd $A/compare
mkdir -p clair3_cp; for f in $A/clair3/*/cp.fasta; do s=$(basename $(dirname $f)); cp $f clair3_cp/$s.fasta; done
for set in bacon_cp_samtools:$A/bacon/cp_samtools/3_assembled/all_assemblies bacon_cp_flye:$A/bacon/cp_flye/3_assembled/all_assemblies \
           ptgaul_cp:$A/compare/ptgaul_oriented clair3_cp:$A/compare/clair3_cp; do
  n=${set%%:*}; d=${set#*:}; [ -d $d ] || continue
  python assembly_score.py $d > scores_assembly_$n.tsv && echo "$n scored"
done
