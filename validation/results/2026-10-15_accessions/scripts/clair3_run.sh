#!/usr/bin/env bash
# Mapping-and-calling pipeline: each accession's R9.4.1 reads mapped to TAIR10.1 with its nuclear copy of the
# mitogenome masked, as for the Illumina truth (reads of other nuclear copies map to the nucleus) with minimap2, Clair3 2.0.3 (model r941_prom_hac_g360+g422,
# --haploid_precise) on the plastome and the mitogenome, then a consensus per organelle: the reference with the
# PASS substitutions, N where the depth is below 10 (any mapping quality, with --min_mq=0: the inverted repeat and the mitogenome repeats have MAPQ 0, as in the Illumina truth). Two accessions at a time, 16 threads each;
# starts after the Illumina truth mapping.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_clair3
B=${BENCH:?the benchmark folder}; W=$B/acc/clair3; cd $W; R=$B/refs/TAIR10.1.orgmask.fasta; M=$CONDA_PREFIX/bin/models/r941_prom_hac_g360+g422
until grep -q DONE $B/acc/truth/illumina.log; do sleep 120; done
one() { f=$1; s=$(basename $f | sed 's/_ERR[0-9]*.fastq.gz//'); [ -s $s/cp.fasta ] && [ -s $s/mt.fasta ] && { echo "$s done"; return; }
  mkdir -p $s
  { [ -s $s/organelles.bam ] || { minimap2 -ax map-ont -t 16 $R $f 2> $s/minimap2.log \
      | samtools view -b -L <(printf 'NC_000932.1\t0\t154478\nNC_037304.1\t0\t367808\n') - | samtools sort -@ 4 -o $s/organelles.bam - \
      && samtools index $s/organelles.bam; }; } \
  && run_clair3.sh --bam_fn=$s/organelles.bam --ref_fn=$R --threads=16 --platform=ont --model_path=$M --output=$s/clair3 \
       --ctg_name=NC_000932.1,NC_037304.1 --haploid_precise --include_all_ctgs --min_mq=0 > $s/clair3.log 2>&1 \
  && for o in cp:NC_000932.1 mt:NC_037304.1; do n=${o%%:*}; c=${o#*:}
       samtools depth -a -Q 0 -r $c $s/organelles.bam | awk '$3 < 10 {print $1"\t"$2-1"\t"$2}' > $s/$n.lowdepth.bed
       bcftools view -f PASS -v snps -r $c -Oz -o $s/$n.snps.vcf.gz $s/clair3/merge_output.vcf.gz && bcftools index -f $s/$n.snps.vcf.gz
       samtools faidx $R $c | bcftools consensus -m $s/$n.lowdepth.bed $s/$n.snps.vcf.gz 2> /dev/null | sed "1s/.*/>$s/" > $s/$n.fasta
     done && echo "$s ok" || echo "$s FAILED"; }
export -f one; export R M
ls $B/accessions/ont/*.fastq.gz | xargs -P 2 -I{} bash -c 'cd '"$W"' && one {}'
echo DONE
