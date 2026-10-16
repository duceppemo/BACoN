#!/usr/bin/env python3
"""Per-position base counts of one sample's Illumina chrM alignment (1000 Genomes 30x CRAM, region chrM), as the
accessions' truth (record 2026-10-15_minimap2, acc_pileup_calls.py): samtools mpileup -B -Q 20 -q 0 -a, written as
<sample>.calls.tsv.gz: position, A, C, G, T counts (deletions and insertions are not counted).

Usage: pileup_calls.py BAM chrM.fasta OUT.tsv.gz"""
import gzip
import re
import subprocess
import sys

INDEL = re.compile(r"[+-](\d+)")


def counts(bases: str) -> dict[str, int]:
    """A, C, G, T in an mpileup bases column ('.'/',' are the reference: 'R')"""
    out = {"A": 0, "C": 0, "G": 0, "T": 0, "R": 0}
    i = 0
    while i < len(bases):
        c = bases[i]
        if c == "^":  # Read start: the next character is its mapping quality
            i += 2
            continue
        if c in "+-":  # Indel: its length, then its bases
            m = INDEL.match(bases, i)
            i = m.end() + int(m.group(1))  # type: ignore[union-attr]
            continue
        u = c.upper()
        if u in out:
            out[u] += 1
        elif c in ".,":
            out["R"] += 1
        i += 1
    return out


def main() -> None:
    bam, reference, out = sys.argv[1:4]
    cmd = ["samtools", "mpileup", "-B", "-Q", "20", "-q", "0", "-a", "-d", "0", "-f", reference, "-r", "chrM", bam]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    with gzip.open(out, "wt") as fh:
        for line in proc.stdout:  # type: ignore[union-attr]
            f = line.split("\t")
            ref = f[2].upper()
            c = counts(f[4])
            if ref in "ACGT":
                c[ref] += c["R"]
            fh.write(f"{f[1]}\t{c['A']}\t{c['C']}\t{c['G']}\t{c['T']}\n")
    if proc.wait():
        sys.exit(f"samtools mpileup failed on {bam}")


if __name__ == "__main__":
    main()
