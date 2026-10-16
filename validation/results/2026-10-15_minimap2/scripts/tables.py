#!/usr/bin/env python3
"""The tables of this record from the runs: the 24 accessions' distances against the Illumina truth (SKA2,
minimap2 and the masking variants tested), the Col-0 lines, and the agreement of templated and Flye assemblies on
the potato data. BENCH: the benchmark folder (acc/, col0/, potato_fref/, potato_mito/); VAL: the folder of
mm2.py's runs (real2/) and of the masking tests (window/). Writes the .tsv files of this record into OUT.

Usage: BENCH=... VAL=... tables.py OUT"""
import os
import sys
from pathlib import Path

BENCH, VAL = Path(os.environ["BENCH"]), Path(os.environ["VAL"])
OUT = Path(sys.argv[1])


def matrix(path: Path) -> tuple[list[str], dict]:
    rows = [line.split("\t") for line in path.read_text().splitlines() if line.strip()]
    names = rows[0][1:]
    return names, {(r[0], b): int(v) for r in rows[1:] for b, v in zip(names, r[1:], strict=True)}


def pairs(names: list[str], reference: bool = False) -> list[tuple[str, str]]:
    names = [n for n in names if reference or n != "Reference"]
    return [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]


def score(path: Path, organelle: str) -> list:
    """Pairs, exact, within 5, mean absolute error, mean signed error, largest error against the truth"""
    t_names, truth = matrix(BENCH / "acc" / "truth" / f"{organelle}.distances.tsv")
    names, m = matrix(path)
    ps = pairs(sorted(set(names) & set(t_names)))
    err = [m[p] - truth[p] for p in ps]
    n = len(err)
    return [n, sum(e == 0 for e in err), sum(abs(e) <= 5 for e in err), round(sum(map(abs, err)) / n, 2),
            round(sum(err) / n, 2), max(map(abs, err))]


def sites(path: Path) -> int:
    lines = path.read_text().split(">")[1:]
    return len("".join(lines[0].splitlines()[1:])) if lines else 0


def write(name: str, header: list[str], rows: list[list]) -> None:
    (OUT / name).write_text("\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))


def main() -> None:
    head = ["pairs", "exact", "within_5", "mean_abs_error", "mean_signed_error", "max_abs_error"]
    rows = []
    for run in ("cp_samtools", "cp_flye", "mt_samtools_fref", "mt_flye_fref"):
        o = run[:2]
        ska = BENCH / "acc" / "bacon" / run / "4_compared" / "ska"
        rows.append([run, "ska", sites(ska / "ska.snps.fasta"), *score(ska / "snp_distances.tsv", o)])
        mm2 = VAL / "real2" / f"acc_{run}"
        rows.append([run, "minimap2", sites(mm2 / "minimap2.snps.fasta"), *score(mm2 / "snp_distances.tsv", o)])
    for o in ("cp", "mt"):
        rows.append([f"{o}_reads", "clair3", "",
                     *score(BENCH / "acc" / "compare" / "distances" / f"clair3_{o}.tsv", o)])
    write("acc_distances.tsv", ["run", "method", "snp_sites", *head], rows)

    # Masking variants (window/window.py, window2.py): w<W>: N within W of a base; w<W>i: N or an indel; d<W>k<K>:
    # at least K N within W; e<W>k<K>: at least K indel events within W (on top of the N rule)
    tests = BENCH / "acc" / "proto" / "window_tests"
    rows = []
    for f in sorted(tests.glob("*.tsv")):
        stem = f.stem  # proto_mm2T_w5_cp, proto_e_F_e50k3_mt
        organelle = stem[-2:]
        assembly = "templated" if "mm2T_" in stem or "_e_T_" in stem else "flye"
        variant = stem.rsplit("_", 2)[-2]
        rows.append([organelle, assembly, variant, *score(f, organelle)])
    write("acc_masking_tests.tsv", ["organelle", "assembly", "variant", *head], rows)

    # Col-0 lines: SNP sites of each method, largest distance between lines and to the reference
    rows = []
    for run in ("cp_samtools", "cp_flye", "cp_samtools_fref", "cp_flye_fref", "mt_samtools", "mt_flye",
                "mt_samtools_fref", "mt_flye_fref"):
        for method, folder in (("ska", BENCH / "col0" / "bacon" / run / "4_compared" / "ska"),
                               ("minimap2", VAL / "real2" / f"col0_{run}")):
            names, m = matrix(folder / "snp_distances.tsv")
            lines = [n for n in names if n != "Reference"]
            snps = folder / ("ska.snps.fasta" if method == "ska" else "minimap2.snps.fasta")
            rows.append([run, method, sites(snps), max(m[p] for p in pairs(names)),
                         max(m["Reference", n] for n in lines)])
    write("col0_lines.tsv", ["run", "method", "snp_sites", "max_between_lines", "max_to_reference"], rows)

    # Potato: templated vs Flye, the same 28 cultivars
    rows = []
    for organelle, base in (("plastome", BENCH / "potato_fref"), ("mitogenome", BENCH / "potato_mito")):
        tag = "potato" if organelle == "plastome" else "potatomt"
        for method, t_path, f_path in (
                ("ska", base / "samtools" / "4_compared" / "ska" / "snp_distances.tsv",
                 base / "flye" / "4_compared" / "ska" / "snp_distances.tsv"),
                ("minimap2", VAL / "real2" / f"{tag}_samtools" / "snp_distances.tsv",
                 VAL / "real2" / f"{tag}_flye" / "snp_distances.tsv")):
            names, t = matrix(t_path)
            _, f = matrix(f_path)
            d = [f[p] - t[p] for p in pairs(names, reference=True)]  # As compare_runs.py: the reference too
            rows.append([organelle, method, len(d), sum(e == 0 for e in d), max(map(abs, d)),
                         round(sum(map(abs, d)) / len(d), 1)])
    write("potato_templated_vs_flye.tsv", ["organelle", "method", "pairs", "identical", "max_difference",
                                           "mean_abs_difference"], rows)


if __name__ == "__main__":
    main()
