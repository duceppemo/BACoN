#!/usr/bin/env bash
# Templated consensus as BACoN makes it (samtools consensus -a -X r10.4_sup --show-ins yes --show-del no -d 3),
# plus -A (IUPAC codes at heterozygous calls), from each sample's kept alignment; ambiguity codes then N, as BACoN
# does (seqio.acgtn).
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate BACoN_eval
A=/media/8tb/analyses/bacon_benchmark/acc; W=$A/proto/ambig
for run in cp_samtools mt_samtools_fref; do mkdir -p $W/$run
  for d in $A/bacon/$run/3_assembled/*/; do s=$(basename $d); [ -f $d/aligned.bam ] || continue
    samtools consensus -a -A -X r10.4_sup --show-ins yes --show-del no -d 3 -f fasta $d/aligned.bam 2>/dev/null \
      | awk -v s=$s 'NR==1{print ">"s; next} /^>/{next} {gsub(/[^ACGTacgt]/,"N"); print toupper($0)}' > $W/$run/$s.fasta
  done; done
