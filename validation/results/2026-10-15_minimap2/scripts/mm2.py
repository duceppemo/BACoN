#!/usr/bin/env python3
"""--snp-method minimap2 (bacon.compare.run_minimap2, core: min_freq 1) on the assemblies of existing BACoN runs,
without touching those runs. Usage: mm2.py OUT RUN_DIR [MIN_FREQ]"""
import sys
import time
from pathlib import Path

from bacon.compare import minimap2_vcf, run_minimap2, snp_distances, write_distances
from bacon.seqio import read_records

out, run = Path(sys.argv[1]), Path(sys.argv[2])
freq = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
assemblies = {p.stem: p for p in sorted((run / "3_assembled" / "all_assemblies").glob("*.fasta"))}
logs = out.parent / (out.name + ".logs")
logs.mkdir(parents=True, exist_ok=True)
t = time.time()
_, snps = run_minimap2(run / "reference.fasta", assemblies, out, logs, threads=16, min_freq=freq)
names, matrix = snp_distances(list(read_records(snps)))
write_distances(out / "snp_distances.tsv", names, matrix)
records = minimap2_vcf(out, run / "reference.fasta", out / "snps.vcf", "test")
sites = len(next(iter(read_records(snps))).seq)
print(f"{out.name}\t{len(assemblies)} genomes\t{sites} sites\t{records} VCF records\t{time.time() - t:.0f} s")
