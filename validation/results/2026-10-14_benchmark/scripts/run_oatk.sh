#!/usr/bin/env bash
# Oatk 1.0 on the six Col-0 lines' R10.4.1 reads, k 501 (the default 1001 gave no mitogenome on MA32_99; see
# run_trial.sh), minimum k-mer coverage 30, embryophyta HMMs of OatkDB v20230921. Three lines at a time.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_truth
B=${BENCH:?the benchmark folder}; W=$B/col0/compare/oatk; mkdir -p $W; cd $W
export DB=$B/tools/OatkDB/v20230921
one() { f=$1; s=$(basename $f | cut -d_ -f1,2); [ -s ${s}_k501.pltd.ctg.fasta ] && { echo "$s done"; return; }
  /usr/bin/time -f "%e s, %M kB" -o ${s}_k501.time oatk -k 501 -c 30 -t 16 -m $DB/embryophyta_mito.fam -p $DB/embryophyta_pltd.fam \
      -o ${s}_k501 $f > ${s}_k501.log 2>&1 && echo "$s ok" || echo "$s FAILED"; }
export -f one
ls $B/col0_lines/ont/*.fastq.gz | xargs -P 3 -I{} bash -c 'cd '"$W"' && one {}'
echo DONE
