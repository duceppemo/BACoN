#!/usr/bin/env bash
# Wall time and memory on one line (MA32_119), each pipeline alone on the machine with 16 threads, one after the
# other: the six-line runs of ptGAUL, TIPPo, CLAW and the Oatk trial overlapped on the machine. Same options as the
# six-line runs. Starts after BACoN's timing runs (run_bacon_timing.sh).
source "$(conda info --base)/etc/profile.d/conda.sh"
B=${BENCH:?the benchmark folder}; C=$B/col0/compare; W=$C/solo; mkdir -p $W; cd $W
until grep -q DONE $C/bacon_timing/run.log; do sleep 120; done
s=MA32_119; f=$(ls $B/col0_lines/ont/${s}_*.fastq.gz); T() { /usr/bin/time -f "%e s, %M kB" -o "$@"; }
go() { echo "[$(date +%T)] $1"; }
go ptgaul; conda activate bench_ptgaul; rm -rf ptgaul; mkdir -p ptgaul/$s
T ptgaul.time ptGAUL.sh -r $B/refs/NC_000932.1.fasta -l $f -t 16 -g 154478 -o ptgaul/$s > ptgaul.log 2>&1 || echo FAILED
go tippo; conda activate bench_tipp; rm -rf tippo; mkdir tippo
(cd tippo && T ../tippo.time TIPPo.v2.4.pl -f $f -p onthq -g organelle -t 16 > ../tippo.log 2>&1) || echo FAILED
go oatk; conda activate bench_truth; D=$B/tools/OatkDB/v20230921
T oatk.time oatk -k 501 -c 30 -t 16 -m $D/embryophyta_mito.fam -p $D/embryophyta_pltd.fam -o oatk_$s $f > oatk.log 2>&1 \
  || echo FAILED
go claw; conda activate bench_claw; rm -rf claw
# A copy of the CLAW clone with this line's reads only; the conda environments of the six-line run
rsync -a --exclude .snakemake --exclude chloro_assembly $C/claw/CLAW/ claw/; mkdir -p claw/chloro_assembly/reads
ln -s $f claw/chloro_assembly/reads/$s.fastq.gz
(cd claw && T ../claw.time snakemake --profile profiles/local --cores 16 --conda-prefix $C/claw/CLAW/.snakemake/conda \
  > ../claw.log 2>&1) || echo FAILED
go DONE
