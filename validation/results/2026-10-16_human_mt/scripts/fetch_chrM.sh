#!/usr/bin/env bash
# Mitochondrial reads of the human benchmark samples, fetched by region from the public alignments:
# ONT R10.4.1 (1000 Genomes long-read consortium; GIAB, ONT open data) and the 1000 Genomes 30x Illumina CRAMs (truth).
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate BACoN_eval
cd "$(dirname "$0")"; mkdir -p ont illumina
# The CRAMs need their reference: sequences fetched by MD5 from EBI, cached here
export REF_PATH="https://www.ebi.ac.uk/ena/cram/md5/%s" REF_CACHE="$PWD/ref_cache/%2s/%2s/%s"
ONT=https://1000g-ont.s3.amazonaws.com/PROCESSED_DATA/ALIGNED_TO_HG38/MINIMAP2_ALIGNED_BAMS
list=$(curl -s "https://1000g-ont.s3.amazonaws.com/?list-type=2&prefix=PROCESSED_DATA/ALIGNED_TO_HG38/MINIMAP2_ALIGNED_BAMS/&max-keys=1000" | grep -o "<Key>[^<]*R10[^<]*\.bam</Key>" | sed 's/<[^>]*>//g;s#.*/##')
while IFS=$'\t' read -r gm na sp pop cram; do
  [ -s ont/$na.fastq.gz ] || { bam=$(echo "$list" | grep "^$gm-"); samtools view -b "$ONT/$bam" chrM 2>>fetch.err | samtools fastq -F 0x900 - 2>>fetch.err | gzip > ont/$na.fastq.gz; }
  [ "$(samtools view -c illumina/$na.chrM.bam 2>/dev/null)" -gt 0 ] 2>/dev/null || samtools view -b -o illumina/$na.chrM.bam "${cram/ftp:\/\//https://}" chrM 2>>fetch.err
  echo "$na $(zcat ont/$na.fastq.gz | awk 'NR%4==2' | wc -l) ONT reads; $(samtools view -c illumina/$na.chrM.bam 2>/dev/null) Illumina reads"
done < samples_1kgp.tsv
GIAB=https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup
for s in HG001:PAW79146 HG002:PAW70337 HG003:PAY87794 HG004:PAY87778 HG005:PAW87816 HG006:PAY77227 HG007:PAY12990; do
  id=${s%%:*}; fc=${s##*:}
  [ -s ont/$id.fastq.gz ] || samtools view -b "$GIAB/$id/$fc/calls.sorted.bam" chrM 2>>fetch.err | samtools fastq -F 0x900 - 2>>fetch.err | gzip > ont/$id.fastq.gz
  echo "$id $(zcat ont/$id.fastq.gz | awk 'NR%4==2' | wc -l) ONT reads"
done
echo DONE
