#!/usr/bin/env python3
"""Mixed-site test, step 2: A, C, G, T counts at every reference position of a sample's templated alignment
(BACoN's aligned.bam), from samtools mpileup -B -Q 0 -q 0 -a (every base: Nanopore base qualities are not
used by the consensus model either).
Usage: mask_ont_counts.py BAM REFERENCE.fasta OUT.tsv.gz"""
import gzip
import subprocess
import sys
from pathlib import Path

# The mpileup parser of the Illumina truth (record 2026-10-15_minimap2)
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "2026-10-15_minimap2" / "scripts"))
from acc_pileup_calls import counts  # noqa: E402

bam, reference, out = sys.argv[1:4]
proc = subprocess.Popen(["samtools", "mpileup", "-B", "-Q", "0", "-q", "0", "-a", "-d", "0", "-f", reference, bam],
                        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
with gzip.open(out, "wt") as fh:
    for line in proc.stdout:
        f = line.split("\t")
        ref = f[2].upper()
        c = counts(f[4])
        if ref in "ACGT":
            c[ref] += c["R"]
        fh.write(f"{f[0]}\t{f[1]}\t{c['A']}\t{c['C']}\t{c['G']}\t{c['T']}\n")
if proc.wait():
    sys.exit(f"samtools mpileup failed on {bam}")
