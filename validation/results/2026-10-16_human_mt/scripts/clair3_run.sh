#!/usr/bin/env bash
# Mapping-and-calling baseline: each sample's chrM reads aligned to rCRS (minimap2 -ax map-ont), Clair3 2.0.3
# --haploid_precise with the R10.4.1 SUP model of its basecaller (v5.0.0, or not stated: r1041_e82_400bps_sup_v500;
# v3.3 and v4.3: r1041_e82_400bps_sup_v410, the closest available), then a consensus: rCRS with the PASS
# substitutions, N where the depth is below 10. Four samples at a time, 8 threads each.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_clair3
B=${BENCH:?the benchmark folder}; H=$B/human_mt; W=$H/clair3; cd $W; R=$H/truth/chrM.fasta; MD=$CONDA_PREFIX/bin/models
[ -f $R.fai ] || samtools faidx $R
one() { s=$1; model=$2; [ -s $s/chrM.fasta ] && { echo "$s done"; return; }
  m=$MD/r1041_e82_400bps_sup_$([ "$model" = v33 ] || [ "$model" = v430 ] && echo v410 || echo v500)
  mkdir -p $s
  { [ -s $s/chrM.bam ] || { minimap2 -ax map-ont -t 8 $R $H/ont/$s.fastq.gz 2> $s/minimap2.log \
      | samtools sort -@ 2 -o $s/chrM.bam - && samtools index $s/chrM.bam; }; } \
  && run_clair3.sh --bam_fn=$s/chrM.bam --ref_fn=$R --threads=8 --platform=ont --model_path=$m --output=$s/clair3 \
       --ctg_name=chrM --haploid_precise --include_all_ctgs --min_mq=0 > $s/clair3.log 2>&1 \
  && samtools depth -a -Q 0 -r chrM $s/chrM.bam | awk '$3 < 10 {print $1"\t"$2-1"\t"$2}' > $s/lowdepth.bed \
  && bcftools view -f PASS -v snps -Oz -o $s/snps.vcf.gz $s/clair3/merge_output.vcf.gz && bcftools index -f $s/snps.vcf.gz \
  && bcftools consensus -f $R -m $s/lowdepth.bed $s/snps.vcf.gz 2> /dev/null | sed "1s/.*/>$s/" > $s/chrM.fasta \
  && echo "$s ok ($(basename $m))" || echo "$s FAILED"; }
export -f one; export H R MD
tail -n +2 $H/samples_ont.tsv | awk -F'\t' '{print $1, $6}' | xargs -P 4 -n 2 bash -c 'cd '"$W"' && one "$0" "$1"'
echo DONE
