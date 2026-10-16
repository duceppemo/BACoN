#!/usr/bin/env python3
"""SKA2 core SNP distances of a folder of genomes and the reference, as BACoN compares genomes
(bacon.compare.run_ska, k 31, sites in every genome). Usage: ska.py REFERENCE.fasta FOLDER OUT.tsv"""
import sys
import tempfile
from pathlib import Path

from bacon.compare import run_ska, snp_distances, write_distances
from bacon.seqio import read_records


def main() -> None:
    reference, folder, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    genomes = {p.stem: p for p in sorted(folder.glob("*.fasta"))}
    with tempfile.TemporaryDirectory(dir=out.parent) as tmp:
        aln, _ = run_ska(reference, genomes, Path(tmp) / "ska", Path(tmp), threads=8, min_freq=1.0)
        names, matrix = snp_distances(list(read_records(aln)))
    write_distances(out, names, matrix)


if __name__ == "__main__":
    main()
