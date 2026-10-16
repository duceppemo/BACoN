#!/usr/bin/env python3
"""Truth of the human samples from their Illumina base counts (pileup_calls.py), as the accessions' (record
2026-10-15_minimap2, acc_truth.py): at each chrM position, a sample's base is the one carried by >= 90% of its
reads at a depth >= 20, else N (low depth, or mixed: heteroplasmy, nuclear copies). Substitutions only.

Writes truth/sites.tsv (every position where any sample has a base other than the reference, or N),
truth/callable.tsv (callable positions and SNPs to the reference per sample) and truth/distances.tsv (pairwise:
positions callable in both samples where their bases differ). Usage: BENCH=... truth.py"""
import gzip
import os
from pathlib import Path

T = Path(os.environ.get("BENCH", ".")) / "human_mt" / "truth"
MIN_DEPTH, MIN_FRACTION = 20, 0.9


def calls(table: Path, length: int) -> list[str]:
    out = ["N"] * length
    with gzip.open(table, "rt") as fh:
        for line in fh:
            pos, *n = line.split("\t")
            depths = [int(x) for x in n]
            total, best = sum(depths), max(range(4), key=depths.__getitem__)
            if total >= MIN_DEPTH and depths[best] >= MIN_FRACTION * total:
                out[int(pos) - 1] = "ACGT"[best]
    return out


def main() -> None:
    ref = "".join(line.strip() for line in (T / "chrM.fasta").read_text().splitlines()[1:]).upper()
    tables = sorted((T / "pileup").glob("*.calls.tsv.gz"))
    samples = [t.name.removesuffix(".calls.tsv.gz") for t in tables]
    bases = {s: calls(t, len(ref)) for s, t in zip(samples, tables, strict=True)}
    variable = []
    with open(T / "sites.tsv", "w") as out:
        out.write("position\treference\t" + "\t".join(samples) + "\n")
        for p in range(len(ref)):
            row = [bases[s][p] for s in samples]
            if any(b != ref[p] for b in row):
                variable.append(row)
                out.write(f"{p + 1}\t{ref[p]}\t" + "\t".join(row) + "\n")
    with open(T / "callable.tsv", "w") as out:
        out.write("sample\tcallable\tnon_reference\n")
        for s in samples:
            snps = sum(b != "N" and b != ref[p] for p, b in enumerate(bases[s]))
            out.write(f"{s}\t{sum(b != 'N' for b in bases[s])}\t{snps}\n")
    with open(T / "distances.tsv", "w") as out:
        out.write("truth\t" + "\t".join(samples) + "\n")
        for i, a in enumerate(samples):
            row = [sum(r[i] != "N" and r[j] != "N" and r[i] != r[j] for r in variable)
                   for j in range(len(samples))]
            out.write(a + "\t" + "\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
