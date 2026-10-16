#!/usr/bin/env python3
"""Reference-coordinate consensus sequences (samtools consensus --show-ins no --show-del yes) of a folder and the
reference as one alignment (pairwise distances: snp-dists) and its core (positions with A, C, G or T in every
sequence). Usage: build_aln.py REFERENCE.fasta FOLDER OUT_PREFIX"""
import sys
from pathlib import Path


def read_one(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines() if not line.startswith(">")).upper()


reference, folder, prefix = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
seqs = {"Reference": read_one(reference)}
seqs.update({f.stem: read_one(f) for f in sorted(folder.glob("*.fasta"))})
core = [i for i in range(len(seqs["Reference"])) if all(v[i] in "ACGT" for v in seqs.values())]
with open(f"{prefix}.aln", "w") as out, open(f"{prefix}.core.aln", "w") as core_out:
    for k, v in seqs.items():
        out.write(f">{k}\n{v}\n")
        core_out.write(f">{k}\n{''.join(v[i] for i in core)}\n")
print(prefix, "core positions", len(core), "of", len(seqs["Reference"]))
