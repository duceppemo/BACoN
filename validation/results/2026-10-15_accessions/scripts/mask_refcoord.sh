#!/usr/bin/env bash
# Mixed-site test, step 1 (reconstructed from the commands run): BACoN's templated consensus in reference
# coordinates (deletions as N, insertions left out: --show-ins no --show-del yes; otherwise BACoN's options) from
# each sample's kept alignment, then one alignment per organelle with the reference (mask_build_aln.py). Step 2:
# mask_ont_counts.py (Nanopore base counts per position); step 3: mask_masked.py (N where the major base is below a
# fraction; core distances).
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate BACoN_eval
A=${BENCH:?the benchmark folder}/acc; P=$A/proto; cd $P
for o in cp:cp_samtools mt:mt_samtools_fref; do org=${o%%:*}; run=${o#*:}; mkdir -p ${org}_refcoord
  for d in $A/bacon/$run/3_assembled/*/; do s=$(basename $d); [ -f $d/aligned.bam ] || continue
    samtools consensus -a -X r10.4_sup --show-ins no --show-del yes -d 3 -f fasta $d/aligned.bam 2>/dev/null \
      | sed "1s/.*/>$s/" | tr '*' 'N' > ${org}_refcoord/$s.fasta; done
  ref=$BENCH/refs/$([ $org = cp ] && echo NC_000932.1 || echo NC_037304.1).fasta
  python mask_build_aln.py $ref ${org}_refcoord ${org}_refcoord
done
