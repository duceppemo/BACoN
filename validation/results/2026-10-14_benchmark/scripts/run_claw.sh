#!/usr/bin/env bash
# CLAW (git clone, Snakemake 7) on the six Col-0 lines: config.yml set to NC_000932.1, --nano-hq (R10.4.1),
# chloroplast_size 154478, cpus 16; the rest as distributed (3,000 reads of >= 5 kb subsampled; Flye and Unicycler).
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_claw
cd "${BENCH:?the benchmark folder}/col0/compare/claw/CLAW"  # git clone of CLAW, config.yml edited as above
/usr/bin/time -f "%e s, %M kB" -o ../claw.time snakemake --profile profiles/local --cores 16 > ../snakemake.log 2>&1 && echo ok || echo FAILED
echo DONE
