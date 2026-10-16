#!/usr/bin/env python3
"""Per-position base counts of one accession's Illumina alignment (samtools mpileup -B -Q 20 -q 0 -a: base quality
>= 20, any mapping quality; reads in the inverted repeat and in the mitogenome's repeats have MAPQ 0, which
bcftools' genotype model all but ignores: it then emitted no allele where every read agreed) written as
<accession>.calls.tsv.gz: organelle, position, A, C, G, T counts (deletions and insertions are not counted).

Usage: pileup_calls.py BAM REFERENCE.fasta OUT.tsv.gz"""
import gzip
import re
import subprocess
import sys

INDEL = re.compile(r"[+-](\d+)")


def counts(bases: str) -> tuple[int, int, int, int]:
    """A, C, G, T in an mpileup bases column ('.'/',' are the reference, handled by the caller as 'R')"""
    out = {"A": 0, "C": 0, "G": 0, "T": 0, "R": 0}
    i = 0
    while i < len(bases):
        c = bases[i]
        if c == "^":  # Read start: the next character is its mapping quality
            i += 2
            continue
        if c in "+-":  # Indel: its length, then its bases
            m = INDEL.match(bases, i)
            i = m.end() + int(m.group(1))
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
    cmd = ["samtools", "mpileup", "-B", "-Q", "20", "-q", "0", "-a", "-d", "0", "-f", reference,
           "-r", "{}", bam]
    with gzip.open(out, "wt") as fh:
        for contig in ("NC_000932.1", "NC_037304.1"):
            args = [a if a != "{}" else contig for a in cmd]
            proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
            for line in proc.stdout:  # type: ignore[union-attr]
                f = line.split("\t")
                ref = f[2].upper()
                c = counts(f[4])
                if ref in "ACGT":
                    c[ref] += c["R"]
                fh.write(f"{contig}\t{f[1]}\t{c['A']}\t{c['C']}\t{c['G']}\t{c['T']}\n")
            if proc.wait():
                sys.exit(f"samtools mpileup failed on {bam}")


if __name__ == "__main__":
    main()
