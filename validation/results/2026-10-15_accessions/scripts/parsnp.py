#!/usr/bin/env python3
"""Parsnp core SNP distances of a folder of genomes and the reference, as BACoN's --snp-method parsnp
(bacon.compare.run_parsnp, -c). Usage: parsnp.py REFERENCE.fasta FOLDER OUT.tsv"""
import sys
from pathlib import Path

from bacon.compare import run_parsnp, snp_distances, write_distances
from bacon.seqio import read_records


def main() -> None:
    reference, folder, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    genomes = {p.stem: p for p in sorted(folder.glob("*.fasta"))}
    work = out.parent / f".{out.stem}"
    work.mkdir(exist_ok=True)
    aln, _ = run_parsnp(reference, genomes, work / "parsnp", work, threads=16)
    names, matrix = snp_distances(list(read_records(aln)))
    write_distances(out, names, matrix)


if __name__ == "__main__":
    main()
