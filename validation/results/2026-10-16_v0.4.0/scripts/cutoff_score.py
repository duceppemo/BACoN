#!/usr/bin/env python3
"""Score the consensus-cutoff variants made by consc.py: distances against the Illumina truths (accessions, human),
the Col-0 lines (truth: 0 between lines), potato (no truth: distinct genomes; two cultivars split at a hairpin).
CUT: the folder of consc.py's outputs (<dataset>/C<cutoff>/); BENCH: the benchmark folder (truths).
Writes a TSV to stdout. Usage: CUT=... BENCH=... cutoff_score.py"""
import os
import statistics as st
from pathlib import Path

C, B = Path(os.environ["CUT"]), Path(os.environ["BENCH"])


def mat(p: Path) -> tuple[list[str], dict]:
    r = [line.split("\t") for line in p.read_text().splitlines()]
    n = r[0][1:]
    return n, {(x[0], b): int(v) for x in r[1:] for b, v in zip(n, x[1:], strict=True)}


def err(dist: Path, truth: Path) -> list:
    n, m = mat(dist)
    tn, t = mat(truth)
    s = sorted(set(n) & set(tn) - {"Reference"})
    e = [m[a, b] - t[a, b] for i, a in enumerate(s) for b in s[i + 1:]]
    return [round(sum(map(abs, e)) / len(e), 2), round(sum(e) / len(e), 2), sum(x == 0 for x in e), len(e)]


def nmed(d: Path) -> float:
    return st.median(int(line.split()[1]) for line in (d / "n_bases.tsv").read_text().splitlines())


def sites(d: Path) -> int:
    t = (d / "mm2" / "minimap2.snps.fasta").read_text().split(">")[1:]
    return len("".join(t[0].splitlines()[1:]))


def main() -> None:
    print("cutoff\tdataset\tN_median\tsnp_sites\tmean_abs_error\tmean_signed_error\texact_pairs\tpairs"
          "\tnote")
    for c in (10, 15, 20, 30):
        truths = (("acc_cp", B / "acc/truth/cp.distances.tsv"), ("acc_mt", B / "acc/truth/mt.distances.tsv"),
                  ("human", B / "human_mt/truth/distances.tsv"))
        for ds, truth in truths:
            d = C / ds / f"C{c}"
            print(c, ds, nmed(d), sites(d), *err(d / "mm2/snp_distances.tsv", truth), "", sep="\t")
        for ds in ("col0_cp", "col0_mt"):
            d = C / ds / f"C{c}"
            _, x = mat(d / "mm2/snp_distances.tsv")
            print(c, ds, nmed(d), sites(d), "", "", "", "", f"max distance {max(x.values())}", sep="\t")
        d = C / "potato" / f"C{c}"
        n, x = mat(d / "mm2/snp_distances.tsv")
        samples, distinct = [a for a in n if a != "Reference"], []
        for a in samples:
            if not any(x[a, k] == 0 for k in distinct):
                distinct.append(a)
        print(c, "potato", nmed(d), sites(d), "", "", "", "",
              f"{len(distinct)} distinct; Bravo {x['Bravo', '12_22_134']}, "
              f"Otrada {x['Otrada', '12_22_134']} SNPs from 12_22_134", sep="\t")


if __name__ == "__main__":
    main()
