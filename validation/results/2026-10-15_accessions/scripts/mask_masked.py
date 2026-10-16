#!/usr/bin/env python3
"""Mixed-site test, step 3: the reference-coordinate templated consensus (mask_refcoord.sh) with N where the
Nanopore reads are mixed: the major base below FRACTION of the A, C, G, T counts (mask_ont_counts.py, in
proto/mask/<run>/). Core distances (positions with a base in every sequence, as --ska-min-freq 1).
Usage: BENCH=... mask_masked.py cp|mt RUN FRACTION OUT.tsv"""
import gzip
import os
import sys
from pathlib import Path

from bacon.compare import snp_distances, write_distances
from bacon.seqio import Record


def read_one(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines() if not line.startswith(">")).upper()


def main() -> None:
    org, run, fraction, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), Path(sys.argv[4])
    bench = Path(os.environ["BENCH"])
    proto = bench / "acc" / "proto"
    name = {"cp": "NC_000932.1", "mt": "NC_037304.1"}[org]
    ref = read_one(bench / "refs" / f"{name}.fasta")
    seqs = {"Reference": ref}
    masked_total = 0
    for f in sorted((proto / f"{org}_refcoord").glob("*.fasta")):
        seq = list(read_one(f))
        with gzip.open(proto / "mask" / run / f"{f.stem}.tsv.gz", "rt") as fh:
            for line in fh:
                _, pos, *n = line.split("\t")
                counts = [int(x) for x in n]
                total = sum(counts)
                if total and max(counts) < fraction * total and seq[int(pos) - 1] in "ACGT":
                    seq[int(pos) - 1] = "N"
                    masked_total += 1
        seqs[f.stem] = "".join(seq)
    core = [i for i in range(len(ref)) if all(v[i] in "ACGT" for v in seqs.values())]
    names, matrix = snp_distances([Record(k, "".join(v[i] for i in core)) for k, v in seqs.items()])
    write_distances(out, names, matrix)
    print(org, run, fraction, "masked per sample", masked_total // (len(seqs) - 1), "core", len(core))


if __name__ == "__main__":
    main()
