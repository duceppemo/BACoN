#!/usr/bin/env bash
# Truth organelle assemblies of the six Col-0 MA32 lines from their PacBio HiFi reads (PRJEB112076), TIPPo v2.4
# reference-free (organelle mode: plastome and mitogenome), into the current folder (BENCH/col0/truth_tippo).
# BENCH: the folder with col0_lines/hifi/. TIPPo's mitochondrial Flye failed for MA32_69 ("No disjointigs"); its
# mitogenome was assembled from TIPPo's mitochondrial reads with:
#   flye --pacbio-hifi MA32_69/*.organelle/*mitochondrial.fasta.filter.fasta -g 400k --asm-coverage 100 \
#        -o MA32_69/mt_flye -t 16
W=$(pwd)
one() { f=$1; s=$(basename $f | cut -d_ -f1,2); mkdir -p $s; cd $s
  [ -e done ] || { /usr/bin/time -f "%e s, %M kB" -o time.txt TIPPo.v2.4.pl -f $f -p hifi -g organelle -t 16 > tippo.log 2>&1 && touch done; }
  echo "$s $( [ -e done ] && echo ok || echo FAILED)"; }
export -f one
ls $BENCH/col0_lines/hifi/*.fastq.gz | xargs -P 3 -I{} bash -c 'cd '"$W"' && one {}'
echo DONE
