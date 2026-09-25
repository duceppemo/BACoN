#!/usr/bin/env python3
"""Generate the BACoN example: a 30 kb circular "organelle" reference and four simulated nanopore samples.

Each sample carries known SNPs relative to the reference, so the expected SNP distances are known:

    sample   SNPs
    alpha    none (identical to the reference)
    beta     6 SNPs
    gamma    beta's 6 SNPs + 4 of its own
    delta    10 SNPs of its own

Reads are 1-20 kb long with about 1.5% errors (substitutions, insertions, deletions), at ~50x depth, mixed with
20% off-target reads from an unrelated sequence, which baiting must remove. Deterministic (fixed seed).
Standard library only.

Usage: python make_example.py OUTPUT_FOLDER
"""

from __future__ import annotations

import gzip
import random
import sys
from pathlib import Path

GENOME_LENGTH = 30_000
DEPTH = 50
ERROR_RATE = 0.015
OFF_TARGET_FRACTION = 0.2
SEED = 20260924

# Planted SNPs: (position, new base index offset); shared ones are applied to several samples.
SNPS = {"alpha": [], "beta": ["b"], "gamma": ["b", "g"], "delta": ["d"]}
SNP_SETS = {"b": 6, "g": 4, "d": 10}

EXPECTED_DISTANCES = {
    ("alpha", "beta"): 6, ("alpha", "gamma"): 10, ("alpha", "delta"): 10,
    ("beta", "gamma"): 4, ("beta", "delta"): 16, ("gamma", "delta"): 20,
}


def random_sequence(rng: random.Random, length: int, gc: float = 0.37) -> str:
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights=weights, k=length))


def mutate(seq: str, positions: list[int], rng: random.Random) -> str:
    s = list(seq)
    for p in positions:
        s[p] = rng.choice([b for b in "ACGT" if b != s[p]])
    return "".join(s)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def add_errors(seq: str, rng: random.Random) -> str:
    out = []
    for base in seq:
        r = rng.random()
        if r < ERROR_RATE * 0.5:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
        elif r < ERROR_RATE * 0.75:
            out.append(base + rng.choice("ACGT"))
        elif r < ERROR_RATE:
            continue
        else:
            out.append(base)
    return "".join(out)


def read_length(rng: random.Random) -> int:
    return max(1000, min(20_000, int(rng.lognormvariate(8.6, 0.5))))


def simulate(genome: str, off_target: str, rng: random.Random, name: str) -> list[str]:
    records, bases, n = [], 0, 0
    circular = genome + genome
    while bases < DEPTH * len(genome):
        n += 1
        length = read_length(rng)
        if rng.random() < OFF_TARGET_FRACTION:
            start = rng.randrange(len(off_target) - length)
            fragment = off_target[start: start + length]
        else:
            start = rng.randrange(len(genome))
            fragment = circular[start: start + min(length, len(genome))]
            bases += len(fragment)
        if rng.random() < 0.5:
            fragment = revcomp(fragment)
        seq = add_errors(fragment, rng)
        qual = "".join(rng.choice("0123456789:;<=>?@ABC") for _ in seq)
        records.append(f"@{name}_read{n} simulated\n{seq}\n+\n{qual}\n")
    return records


def main(out: Path) -> None:
    rng = random.Random(SEED)
    reference = random_sequence(rng, GENOME_LENGTH)
    off_target = random_sequence(rng, 200_000, gc=0.5)
    positions = rng.sample(range(500, GENOME_LENGTH - 500), sum(SNP_SETS.values()))
    sets, i = {}, 0
    for key, count in SNP_SETS.items():
        sets[key] = sorted(positions[i: i + count])
        i += count
    (out / "reads").mkdir(parents=True, exist_ok=True)
    with open(out / "reference.fasta", "w") as fh:
        fh.write(">organelle simulated circular reference\n")
        for j in range(0, len(reference), 80):
            fh.write(reference[j: j + 80] + "\n")
    with open(out / "planted_snps.tsv", "w") as fh:
        fh.write("sample\tpositions (1-based)\n")
        for sample, keys in SNPS.items():
            fh.write(f"{sample}\t{','.join(str(p + 1) for k in keys for p in sets[k])}\n")
    for sample, keys in SNPS.items():
        genome = reference
        for k in keys:
            genome = mutate(genome, sets[k], random.Random(f"{SEED}-{k}"))
        with gzip.GzipFile(out / "reads" / f"{sample}.fastq.gz", "wb", compresslevel=6, mtime=0) as fh:
            fh.write("".join(simulate(genome, off_target, rng, sample)).encode())
    print(f"Example written to {out}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(Path(sys.argv[1]))
