#!/usr/bin/env python3
"""Distances between the Col-0 lines for each tool's assemblies: SKA2 core SNPs (k 31, as BACoN's default
comparison: bacon.compare.run_ska), with the reference. Writes distances/<run>.tsv and prints a summary line.
Usage: distances.py RUN_FOLDER... (each with 3_assembled/all_assemblies/<line>.fasta)."""
import os
import sys
import tempfile
from pathlib import Path

from bacon.compare import run_ska, snp_distances, write_distances
from bacon.seqio import read_records

REFS = Path(os.environ.get("BENCH", ".")) / "refs"  # BENCH: the folder with refs/
OUT = Path(os.environ.get("BENCH", ".")) / "col0" / "compare" / "distances"


def main() -> None:
    OUT.mkdir(exist_ok=True)
    for run in map(Path, sys.argv[1:]):
        name = run.name if run.name[:3] not in ("cp_", "mt_") else "bacon_" + run.name
        organelle = "cp" if "cp" in name else "mt"
        reference = REFS / ("NC_000932.1.fasta" if organelle == "cp" else "NC_037304.1.fasta")
        assemblies = {p.stem: p for p in sorted((run / "3_assembled" / "all_assemblies").glob("*.fasta"))}
        with tempfile.TemporaryDirectory(dir=OUT) as tmp:
            aln, _ = run_ska(reference, assemblies, Path(tmp) / "ska", Path(tmp), threads=8, min_freq=1.0)
            names, matrix = snp_distances(list(read_records(aln)))
            sites = len(next(iter(read_records(aln))).seq) if assemblies else 0
        write_distances(OUT / f"{name}.tsv", names, matrix)
        i = names.index("Reference")
        to_ref = [matrix[i][j] for j in range(len(names)) if j != i]
        between = [matrix[a][b] for a in range(len(names)) for b in range(a + 1, len(names)) if i not in (a, b)]
        print(f"{name}\tlines {len(assemblies)}\tcore SNP sites {sites}\tto reference {min(to_ref)}-{max(to_ref)}"
              f"\tbetween lines max {max(between) if between else 'NA'}")


if __name__ == "__main__":
    main()
