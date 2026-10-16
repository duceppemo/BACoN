#!/usr/bin/env bash
# Re-run the minimap2 comparisons with the final code (N-density masking); the earlier ones are moved aside.
set -u
V=/media/8tb/analyses/bacon_validation_mm2
for grp in all new; do
  for sc in $(ls -d $V/$grp/*/ | xargs -n1 basename | grep -v sim); do
    for run in samtools flye; do
      d=$V/$grp/$sc/$run
      [ -d $d/4_compared/minimap2 ] || continue
      mv $d/4_compared/minimap2 $d/minimap2_nomask
      extra=$([ $run = flye ] && echo "-a flye")
      echo "[$(date +%T)] $grp/$sc/$run"
      bacon -r $V/$grp/sim/$sc/reference.fasta -i $V/$grp/sim/$sc/reads -o $d -t 16 -p 4 $extra --snp-method minimap2 >> $d.rerun.log 2>&1 || echo "  FAILED"
    done
    python /home/bioinfo/prog/BACoN/validation/score.py $V/$grp/sim/$sc $V/$grp/$sc/*/ --markdown $V/$grp/SCORES_$sc.md --json $V/$grp/SCORES_$sc.json > /dev/null
  done
done
echo DONE
