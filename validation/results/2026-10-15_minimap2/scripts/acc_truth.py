#!/usr/bin/env python3
"""Truth of the 24 accessions from their Illumina reads (acc_truth_illumina.sh: the alignments; its VCFs are not
used): at each organelle position, a sample's base is the one carried by >= 90% of its reads at a depth >= 20
(base counts at any mapping quality, base quality >= 20: acc_pileup_calls.py, into truth/pileup/), else N (low
depth, or mixed: nuclear copies that still map, heteroplasmy, repeat copies that differ). Substitutions only
(SKA2 and the other pipelines' distances count substitutions).

Writes truth/<organelle>.sites.tsv (every position where any sample has a base other than the reference, or N),
truth/<organelle>.callable.tsv (callable positions per sample) and truth/<organelle>.distances.tsv (pairwise:
positions callable in both samples where their bases differ).

Usage: BENCH=... truth.py"""
import gzip
import os
from pathlib import Path

B = Path(os.environ.get("BENCH", ".")) / "acc" / "truth"
REFS = Path(os.environ.get("BENCH", ".")) / "refs"
ORGANELLES = {"NC_000932.1": "cp", "NC_037304.1": "mt"}
MIN_DEPTH, MIN_FRACTION = 20, 0.9


def calls(table: Path) -> dict:
    """{organelle: {position: base or 'N'}} from the base counts of pileup_calls.py"""
    out = {o: {} for o in ORGANELLES.values()}
    with gzip.open(table, "rt") as fh:
        for line in fh:
            contig, pos, *n = line.split("\t")
            depths = [int(x) for x in n]
            total, best = sum(depths), max(range(4), key=depths.__getitem__)
            ok = total >= MIN_DEPTH and depths[best] >= MIN_FRACTION * total
            out[ORGANELLES[contig]][int(pos)] = "ACGT"[best] if ok else "N"
    return out


def main() -> None:
    vcfs = sorted((B / "pileup").glob("*.calls.tsv.gz"))  # pileup_calls.py
    samples = [v.name.removesuffix(".calls.tsv.gz") for v in vcfs]
    lengths = {"cp": 154478, "mt": 367808}
    per_sample = {s: calls(v) for s, v in zip(samples, vcfs)}
    for organelle, length in lengths.items():
        # The reference's base at each position (1-based)
        name = next(k for k, v in ORGANELLES.items() if v == organelle)
        seq = (REFS / f"{name}.fasta").read_text().split("\n", 1)[1].replace("\n", "").upper()
        ref = {i + 1: b for i, b in enumerate(seq)}
        bases = {s: [per_sample[s][organelle].get(p, "N") for p in range(1, length + 1)] for s in samples}
        variable = []  # Rows where any sample differs from the reference or is N: elsewhere all samples agree
        with open(B / f"{organelle}.sites.tsv", "w") as out:
            out.write("position\treference\t" + "\t".join(samples) + "\n")
            for p in range(1, length + 1):
                row = [bases[s][p - 1] for s in samples]
                if any(b != ref.get(p) for b in row):
                    variable.append(row)
                    out.write(f"{p}\t{ref.get(p, 'N')}\t" + "\t".join(row) + "\n")
        with open(B / f"{organelle}.callable.tsv", "w") as out:
            out.write("sample\tcallable\tnon_reference\n")
            for s in samples:
                callable_ = sum(b != "N" for b in bases[s])
                snps = sum(b != "N" and b != ref.get(p + 1) for p, b in enumerate(bases[s]))
                out.write(f"{s}\t{callable_}\t{snps}\n")
        with open(B / f"{organelle}.distances.tsv", "w") as out:
            out.write("truth\t" + "\t".join(samples) + "\n")
            for i, a in enumerate(samples):
                row = [sum(r[i] != "N" and r[j] != "N" and r[i] != r[j] for r in variable)
                       for j in range(len(samples))]
                out.write(a + "\t" + "\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
