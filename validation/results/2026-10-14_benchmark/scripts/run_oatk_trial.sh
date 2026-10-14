#!/usr/bin/env bash
# Oatk 1.0 trial on one Col-0 line's R10.4.1 reads (Oatk is designed for HiFi): the default k-mer (1001) and a
# smaller one (501), minimum k-mer coverage 30, embryophyta HMMs of OatkDB v20230921.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_truth
B=${BENCH:?the benchmark folder}; W=$B/col0/compare/oatk; mkdir -p $W; cd $W
DB=$B/tools/OatkDB/v20230921; f=$(ls $B/col0_lines/ont/MA32_99_*.fastq.gz)
for k in 1001 501; do
  /usr/bin/time -f "%e s, %M kB" -o MA32_99_k$k.time oatk -k $k -c 30 -t 16 -m $DB/embryophyta_mito.fam -p $DB/embryophyta_pltd.fam \
      -o MA32_99_k$k $f > MA32_99_k$k.log 2>&1 && echo "k$k ok" || echo "k$k FAILED"
done
echo DONE
