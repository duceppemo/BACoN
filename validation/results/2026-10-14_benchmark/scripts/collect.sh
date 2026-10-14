#!/usr/bin/env bash
# Gather each pipeline's final assemblies as compare/collect/<tool>_<cp|mt>/3_assembled/all_assemblies/<line>.fasta
# (the layout of a BACoN run, which score.py and distances.py read).
B=${BENCH:?the benchmark folder}; cd $B/col0/compare
put() { mkdir -p collect/$1/3_assembled/all_assemblies; cp "$2" collect/$1/3_assembled/all_assemblies/$3.fasta; }
for s in MA32_119 MA32_59 MA32_69 MA32_79 MA32_89 MA32_99; do
  # ptGAUL: its final assembly (path1: the first of the two orientations of the single-copy region)
  put ptgaul_cp ptgaul/$s/result_3000/ptGAUL_final_assembly/path1.fasta $s
  # TIPPo: the final plastome, or the last round's when the last round already gave one contig; the mitogenome
  # from its Flye run (absent where Flye failed)
  d=$(ls -d tippo_ont/$s/*.organelle); cp=$(ls -t $d/*.edge_*.organelle.chloroplast.fasta 2>/dev/null | head -1)
  [ -z "$cp" ] && cp=$(ls -t $d/*.round*.organelle.chloroplast.fasta | head -1); put tippo_cp $cp $s
  mt=$d/$(basename $d .organelle).mitochondrial.fasta.filter.fasta.flye/assembly.fasta
  [ -s $mt ] && put tippo_mt $mt $s || echo "$s: no TIPPo mitogenome"
  # Oatk, k 501: plastid and mitochondrial contigs
  put oatk_cp oatk/${s}_k501.pltd.ctg.fasta $s; put oatk_mt oatk/${s}_k501.mito.ctg.fasta $s
  # CLAW: the Flye and the Unicycler plastome
  for a in flye unicycler; do f=claw/CLAW/chloro_assembly/$s~${a}_chloroplast.fasta; [ -s $f ] && put claw_${a}_cp $f $s; done
done
