#!/usr/bin/env bash
# Illumina truth of the 24 accessions: the first 10 million read pairs (organelles at hundreds to thousands x)
# mapped to TAIR10.1 (nuclear genome, NC_037304.1, NC_000932.1; reads of nuclear copies of organelle DNA map to
# the nucleus) with its 259 kb copy of the mitogenome on chromosome 2 masked (refs/nuc_mask.bed: 95-99% identical
# to the mitogenome, its reads could not be assigned), then every organelle position piled up at any mapping
# quality (reads in the inverted repeat and in the mitogenome's repeats have MAPQ 0; base quality >= 20) and called
# haploid; truth.py keeps the sites with one clear base (reads of less similar nuclear copies make sites mixed). Four
# accessions at a time, 16 threads each.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_snp
B=${BENCH:?the benchmark folder}; W=$B/acc/truth/illumina; mkdir -p $W; cd $W; R=$B/refs/TAIR10.1.orgmask.fasta
one() { r1=$1; s=$(basename $r1 | sed 's/_ERR[0-9]*_1.fastq.gz//'); r2=${r1/_1.fastq.gz/_2.fastq.gz}
  [ -s $s.vcf.gz ] && { echo "$s done"; return; }
  bwa mem -t 16 -R "@RG\tID:$s\tSM:$s" $R <(zcat $r1 | head -n 40000000) <(zcat $r2 | head -n 40000000) 2> $s.bwa.log | samtools sort -@ 4 -m 2G -o $s.bam - && samtools index $s.bam \
  && bcftools mpileup -f $R -r NC_000932.1,NC_037304.1 -q 0 -Q 20 -d 5000 -a AD,DP -Ou $s.bam 2> $s.mpileup.log \
     | bcftools call -m --ploidy 1 -Oz -o $s.vcf.gz && echo "$s ok" || echo "$s FAILED"; }
export -f one; export R
ls $B/accessions/illumina/*_1.fastq.gz | xargs -P 4 -I{} bash -c 'cd '"$W"' && one {}'
echo DONE
