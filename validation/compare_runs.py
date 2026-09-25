#!/usr/bin/env python3
"""Compare BACoN runs on the same real dataset, where no truth is known.

Usage: python compare_runs.py RUN_FOLDER [RUN_FOLDER ...]

Prints, for each run: samples assembled, contigs, circular contigs, assembly length relative to the reference;
then, for every pair of SNP distance matrices (4_compared/*/snp_distances.tsv of each run): the samples in
common, how many pairwise distances are identical, the largest difference, and the Pearson correlation.
Agreement between independent methods (a de novo and a templated assembly; split k-mers and alignment) is the
best available evidence on real data.
"""

from __future__ import annotations

import statistics
import sys
from itertools import combinations
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line.strip()]


def matrix(path: Path) -> dict[tuple[str, str], int]:
    rows = read_tsv(path)
    first = next(iter(rows[0]))
    out = {}
    for r in rows:
        for k, v in r.items():
            if k != first:
                out[tuple(sorted((r[first], k)))] = int(v)
    return {k: v for k, v in out.items() if k[0] != k[1]}


def pearson(x: list[float], y: list[float]) -> float:
    if len(x) < 3 or statistics.pstdev(x) == 0 or statistics.pstdev(y) == 0:
        return float("nan")
    return statistics.correlation(x, y)


def main(runs: list[Path]) -> None:
    print(f"{'run':22}{'assembled':>10}{'contigs/sample':>16}{'circular':>9}{'len/ref (min-max)':>20}")
    matrices = {}
    for run in runs:
        rows = read_tsv(run / "summary.tsv")
        ok = [r for r in rows if r["Status"] == "ok"]
        contigs = [int(r["Contigs"]) for r in ok if r["Contigs"].isdigit()]
        circ = sum(int(r["Circular_contigs"]) for r in ok if r["Circular_contigs"].isdigit())
        ratios = [float(r["Length_vs_reference"]) for r in ok if r["Length_vs_reference"] not in ("NA", "")]
        span = f"{min(ratios):.3f}-{max(ratios):.3f}" if ratios else "NA"
        mean_contigs = f"{statistics.mean(contigs):.2f}" if contigs else "NA"
        print(f"{run.name:22}{len(ok):>6}/{len(rows):<3}{mean_contigs:>16}{circ:>9}{span:>20}")
        for dist in sorted((run / "4_compared").glob("*/snp_distances.tsv")):
            matrices[f"{run.name}/{dist.parent.name}"] = matrix(dist)
    print()
    print(f"{'A':28}{'B':28}{'pairs':>7}{'identical':>10}{'max diff':>9}{'r':>7}")
    for a, b in combinations(matrices, 2):
        common = sorted(set(matrices[a]) & set(matrices[b]))
        if not common:
            continue
        x = [matrices[a][k] for k in common]
        y = [matrices[b][k] for k in common]
        same = sum(1 for i, j in zip(x, y) if i == j)
        print(f"{a:28}{b:28}{len(common):>7}{same:>10}{max(abs(i - j) for i, j in zip(x, y)):>9}"
              f"{pearson(x, y):>7.3f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main([Path(p) for p in sys.argv[1:]])
