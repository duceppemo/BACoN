#!/usr/bin/env python3
"""Score each pipeline's distances and SNPs against the Illumina truth (truth/truth.py).

Distances: for every pair of accessions present in both matrices, the method's distance minus the truth's
(positions callable in both by Illumina): pairs exact, within 1, within 5, mean absolute error, mean signed error
(negative: distances too low), largest error. Writes compare/scores_distances.tsv and compare/errors/<method>.tsv.

Distances at judgeable positions (methods whose genomes are in reference coordinates: an alignment whose every
column is a reference position): pairwise differences between A, C, G, T at positions callable by Illumina in both
accessions only, so that calls where the truth cannot judge (mixed or low-depth Illumina sites: nuclear copies,
repeats) do not count. Writes compare/scores_distances_judgeable.tsv.

SNPs (methods with SNP calls per sample, in reference coordinates): at the positions callable by Illumina in that
accession, true SNPs found (recall) and calls that are true (precision). Writes compare/scores_snps.tsv.

Usage: BENCH=... score.py"""
import csv
import gzip
import os
from pathlib import Path

A = Path(os.environ.get("BENCH", ".")) / "acc"
C = A / "compare"
NAMES = {"NC_000932.1": "cp", "NC_037304.1": "mt"}


def matrix(path: Path) -> dict:
    rows = [line.rstrip("\n").split("\t") for line in path.read_text().splitlines() if line.strip()]
    names = rows[0][1:]
    return {(r[0], b): int(v) for r in rows[1:] for b, v in zip(names, r[1:])}


def methods() -> dict:
    """{name: (organelle, distance matrix path)}"""
    out = {}
    for run in ("cp_samtools", "cp_flye", "mt_samtools_fref", "mt_flye_fref"):
        out[f"bacon_ska_{run}"] = (run[:2], A / "bacon" / run / "4_compared" / "ska" / "snp_distances.tsv")
    for f in sorted((C / "distances").glob("*.tsv")):
        out[f.stem] = ("cp" if f.stem.endswith("_cp") else "mt", f)
    return {k: v for k, v in out.items() if v[1].is_file()}


def truth_sites(organelle: str) -> tuple[list[str], dict]:
    """Samples and {position: (reference, {sample: base})} of truth/<organelle>.sites.tsv"""
    with open(A / "truth" / f"{organelle}.sites.tsv") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[2:]
        sites = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sites[int(f[0])] = (f[1], dict(zip(samples, f[2:])))
    return samples, sites


def vcf_snps(path: Path, sample_column: str | None = None) -> dict:
    """{sample: {(organelle, position): alt}} from a VCF (one or several samples; haploid or diploid GT)"""
    out = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        samples = []
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:] if sample_column is None else [sample_column]
                continue
            organelle = NAMES.get(f[0])
            ref, alts = f[3].upper(), f[4].upper().split(",")
            if organelle is None or len(ref) != 1 or (len(f) > 6 and f[6] not in ("PASS", ".")):
                continue
            keys = f[8].split(":")
            for s, values in zip(samples, f[9:]):
                gt = dict(zip(keys, values.split(":"))).get("GT", ".").replace("|", "/").split("/")[0]
                if gt.isdigit() and int(gt) > 0 and len(alts[int(gt) - 1]) == 1:
                    out.setdefault(s, {})[(organelle, int(f[1]))] = alts[int(gt) - 1]
    return out


def snp_calls() -> dict:
    """{method: (organelle, {sample: {position: alt}})}"""
    out = {}
    for run in ("cp_samtools", "cp_flye", "mt_samtools_fref", "mt_flye_fref"):
        vcf = A / "bacon" / run / "4_compared" / "ska" / "snps.vcf"
        if vcf.is_file():
            calls = vcf_snps(vcf)
            out[f"bacon_ska_{run}"] = (run[:2], {s: {p: a for (o, p), a in c.items()} for s, c in calls.items()})
    for run in ("cp_samtools", "cp_flye", "mt_samtools_fref", "mt_flye_fref"):
        vcf = C / "minimap2" / run / "snps.vcf"  # --snp-method minimap2: every position where genomes differ
        if vcf.is_file():
            calls = vcf_snps(vcf)
            out[f"bacon_minimap2_{run}"] = (run[:2],
                                            {s: {p: a for (o, p), a in c.items()} for s, c in calls.items()})
    for organelle in ("cp", "mt"):
        clair = {}
        for vcf in sorted((A / "clair3").glob(f"*/{organelle}.snps.vcf.gz")):
            calls = vcf_snps(vcf, vcf.parent.name).get(vcf.parent.name, {})
            clair[vcf.parent.name] = {p: a for (o, p), a in calls.items()}
        if clair:
            out[f"clair3_{organelle}"] = (organelle, clair)
        snippy = {}
        for vcf in sorted((C / f"snippy_{organelle}").glob("*/snps.vcf")):
            calls = vcf_snps(vcf, vcf.parent.name).get(vcf.parent.name, {})
            snippy[vcf.parent.name] = {p: a for (o, p), a in calls.items()}
        if snippy:
            out[f"snippy_{organelle}"] = (organelle, snippy)
    return out


def alignments() -> dict:
    """{name: (organelle, alignment, positions)}: alignments whose columns are reference positions, every position
    (positions None) or those listed (minimap2: every position where genomes differ, minimap2.positions.tsv)"""
    out = {}
    for o in ("cp", "mt"):
        out[f"clair3_{o}"] = (o, C / f"clair3_{o}.aln", None)
        out[f"snippy_{o}"] = (o, C / f"snippy_{o}" / "core.full.aln", None)
    for run in ("cp_samtools", "cp_flye", "mt_samtools_fref", "mt_flye_fref"):
        d = C / "minimap2" / run
        if (d / "minimap2.positions.tsv").is_file():
            lines = (d / "minimap2.positions.tsv").read_text().splitlines()
            positions = [int(line.split("\t")[1]) for line in lines]
            out[f"bacon_minimap2_{run}"] = (run[:2], d / "minimap2.variable.fasta", positions)
    return {k: v for k, v in out.items() if v[1].is_file()}


def read_aln(path: Path) -> dict:
    seqs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def judgeable(organelle: str, aln: Path, truth: dict, sites: tuple, positions: list | None = None) -> list:
    """Pairs (a, b, truth, method) counting differences at positions callable by Illumina in a and b (column k is
    position k + 1, or positions[k])"""
    samples, table = sites
    seqs = read_aln(aln)
    names = sorted(set(seqs) & set(samples))
    # Positions not in the sites table are callable in every accession (all reference): differences there count
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            sa, sb = seqs[a], seqs[b]
            masked = {p for p, v in table.items() if v[1][a] == "N" or v[1][b] == "N"}
            d = sum(1 for k, (x, y) in enumerate(zip(sa, sb, strict=False))
                    if x != y and x in "ACGT" and y in "ACGT"
                    and (positions[k] if positions else k + 1) not in masked)
            rows.append((a, b, truth[a, b], d))
    return rows


def main() -> None:
    truths = {o: matrix(A / "truth" / f"{o}.distances.tsv") for o in ("cp", "mt")}
    (C / "errors").mkdir(exist_ok=True)
    with open(C / "scores_distances.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["method", "organelle", "accessions", "pairs", "exact", "within_1", "within_5",
                    "mean_abs_error", "mean_signed_error", "max_abs_error"])
        for name, (organelle, path) in methods().items():
            m, t = matrix(path), truths[organelle]
            samples = sorted({a for a, _ in m} & {a for a, _ in t})
            pairs = [(a, b) for i, a in enumerate(samples) for b in samples[i + 1:]]
            err = [m[p] - t[p] for p in pairs]
            with open(C / "errors" / f"{name}.tsv", "w") as out:
                out.write("a\tb\ttruth\tmethod\terror\n")
                for (a, b), e in zip(pairs, err):
                    out.write(f"{a}\t{b}\t{t[a, b]}\t{m[a, b]}\t{e}\n")
            n = len(err) or 1
            w.writerow([name, organelle, len(samples), len(pairs), sum(e == 0 for e in err),
                        sum(abs(e) <= 1 for e in err), sum(abs(e) <= 5 for e in err),
                        round(sum(map(abs, err)) / n, 2), round(sum(err) / n, 2), max(map(abs, err), default=0)])
    sites = {o: truth_sites(o) for o in ("cp", "mt")}
    with open(C / "scores_distances_judgeable.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["method", "organelle", "accessions", "pairs", "exact", "within_1", "within_5",
                    "mean_abs_error", "mean_signed_error", "max_abs_error"])
        for name, (organelle, aln, positions) in alignments().items():
            rows = judgeable(organelle, aln, truths[organelle], sites[organelle], positions)
            err = [m - t for _, _, t, m in rows]
            n = len(err) or 1
            w.writerow([name, organelle, len({r[0] for r in rows} | {r[1] for r in rows}), len(rows),
                        sum(e == 0 for e in err), sum(abs(e) <= 1 for e in err), sum(abs(e) <= 5 for e in err),
                        round(sum(map(abs, err)) / n, 2), round(sum(err) / n, 2), max(map(abs, err), default=0)])
    with open(C / "scores_snps.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["method", "organelle", "accessions", "true_snps", "found", "recall", "calls_at_callable",
                    "true_calls", "precision"])
        for name, (organelle, calls) in snp_calls().items():
            samples, table = sites[organelle]
            found = true_snps = called = true_calls = 0
            accessions = [s for s in calls if s in samples]
            for s in accessions:
                truth = {p: v[1][s] for p, v in table.items() if v[1][s] not in ("N", v[0])}
                uncallable = {p for p, v in table.items() if v[1][s] == "N"}
                c = calls[s]
                true_snps += len(truth)
                found += sum(c.get(p) == b for p, b in truth.items())
                at_callable = {p: b for p, b in c.items() if p not in uncallable}
                called += len(at_callable)
                true_calls += sum(truth.get(p) == b for p, b in at_callable.items())
            w.writerow([name, organelle, len(accessions), true_snps, found, round(found / (true_snps or 1), 3),
                        called, true_calls, round(true_calls / (called or 1), 3)])


if __name__ == "__main__":
    main()
