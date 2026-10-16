#!/usr/bin/env bash
# Comparisons of assemblies (after ptgaul/run.sh and bacon/run.sh):
# - plastome: ptGAUL's plastomes, oriented and rotated to the reference (orient.py), then (1) MAFFT --auto of the 24
#   and the reference, snp-dists; (2) Snippy --ctgs on each, snippy-core, snp-dists of core.aln; (3) SKA2 as BACoN
#   compares genomes (ska.py)
# - mitogenome: BACoN's Flye assemblies (--filter-by-reference), Snippy --ctgs, snippy-core, snp-dists
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_snp
B=${BENCH:?the benchmark folder}; W=$B/acc/compare; cd $W; D=$W/distances
until grep -q DONE $B/acc/ptgaul/run.log && grep -q DONE $B/acc/bacon/run.log; do sleep 120; done
mkdir -p ptgaul_oriented
# ptGAUL's paths; where its graph had a single edge it names final_assembly.fasta as the plastome but writes none
# (ptGAUL 1.0.5): then its Flye assembly, which that file would hold
for d in $B/acc/ptgaul/*/; do s=$(basename $d); p=$(ls $d/result_*/ptGAUL_final_assembly/*.fasta 2>/dev/null)
  [ -z "$p" ] && p=$(ls $d/result_*/flye_cpONT/assembly.fasta 2>/dev/null)
  [ -n "$p" ] && python orient.py $B/refs/NC_000932.1.fasta ptgaul_oriented/$s.fasta $s $p || echo "$s: no ptGAUL plastome"; done
# (1) MAFFT
cat $B/refs/NC_000932.1.fasta ptgaul_oriented/*.fasta | sed '1s/.*/>Reference/' > mafft_cp.input.fasta
mafft --auto --thread 16 mafft_cp.input.fasta > mafft_cp.aln 2> mafft_cp.log && snp-dists -q mafft_cp.aln > $D/ptgaul_mafft_cp.tsv
# (2) Snippy, plastome and mitogenome
snip() { org=$1; ref=$2; shift 2; mkdir -p snippy_$org
  for f in "$@"; do s=$(basename $f .fasta); [ -s snippy_$org/$s/snps.aligned.fa ] || \
    snippy --ctgs $f --ref $ref --outdir snippy_$org/$s --cpus 16 --force > snippy_$org/$s.log 2>&1 || echo "snippy $org $s FAILED"; done
  (cd snippy_$org && snippy-core --ref $ref --prefix core $(ls -d */) > core.log 2>&1) && snp-dists -q snippy_$org/core.aln > $D/snippy_$org.tsv; }
snip cp $B/refs/NC_000932.1.fasta ptgaul_oriented/*.fasta
snip mt $B/refs/NC_037304.1.fasta $B/acc/bacon/mt_flye_fref/3_assembled/all_assemblies/*.fasta
# (3) SKA2 on the ptGAUL plastomes
conda activate BACoN_eval; python ska.py $B/refs/NC_000932.1.fasta ptgaul_oriented $D/ptgaul_ska_cp.tsv
echo DONE
