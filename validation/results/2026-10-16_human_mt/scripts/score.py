#!/usr/bin/env python3
"""Score each pipeline on the human samples against the Illumina truth (truth/truth.py), as the accessions'
score.py (record 2026-10-15_accessions):

- distances: every pair present in both matrices, the method's distance minus the truth's (exact, within 1, within
  5, mean absolute error, mean signed error, largest error), and the pairs that must be 0 (mother and child;
  identical in the truth). Writes scores_distances.tsv, errors/<method>.tsv;
- distances at judgeable positions (sequences in reference coordinates): differences only at positions the truth
  calls in both samples. Writes scores_distances_judgeable.tsv;
- SNPs per sample (VCFs): recall and precision at positions the truth calls in that sample. Writes scores_snps.tsv.

Usage: BENCH=... score.py"""
import csv
import gzip
import os
from pathlib import Path

H = Path(os.environ.get("BENCH", ".")) / "human_mt"
RUNS = ("samtools", "flye")
ZERO = [("HG02024", "HG02025")]  # Mother and child: the same mtDNA


def matrix(path: Path) -> dict:
    rows = [line.split("\t") for line in path.read_text().splitlines() if line.strip()]
    names = rows[0][1:]
    return {(r[0], b): int(v) for r in rows[1:] for b, v in zip(names, r[1:], strict=True)}


def methods() -> dict:
    out = {}
    for run in RUNS:
        out[f"bacon_ska_{run}"] = H / "bacon" / run / "4_compared" / "ska" / "snp_distances.tsv"
        out[f"bacon_minimap2_{run}"] = H / "bacon" / run / "4_compared" / "minimap2" / "snp_distances.tsv"
    out["clair3"] = H / "compare" / "clair3.distances.tsv"
    return {k: v for k, v in out.items() if v.is_file()}


def truth_sites() -> tuple[list[str], dict]:
    with open(H / "truth" / "sites.tsv") as fh:
        samples = fh.readline().rstrip("\n").split("\t")[2:]
        sites = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sites[int(f[0])] = (f[1], dict(zip(samples, f[2:], strict=True)))
    return samples, sites


def vcf_snps(path: Path, sample: str | None = None) -> dict:
    """{sample: {position: alt}} (single-base REF and ALT; PASS or '.')"""
    out: dict = {}
    with (gzip.open if path.suffix == ".gz" else open)(path, "rt") as fh:
        samples: list[str] = []
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:] if sample is None else [sample]
                continue
            alts = f[4].upper().split(",")
            if len(f[3]) != 1 or f[6] not in ("PASS", "."):
                continue
            keys = f[8].split(":")
            for s, values in zip(samples, f[9:], strict=False):
                fields = dict(zip(keys, values.split(":"), strict=False))
                gt = fields.get("GT", ".").replace("|", "/").split("/")[0]
                if gt.isdigit() and int(gt) > 0 and len(alts[int(gt) - 1]) == 1:
                    out.setdefault(s, {})[int(f[1])] = alts[int(gt) - 1]
    return out


def snp_calls() -> dict:
    out = {}
    for run in RUNS:
        for method in ("ska", "minimap2"):
            vcf = H / "bacon" / run / "4_compared" / method / "snps.vcf"
            if vcf.is_file():
                out[f"bacon_{method}_{run}"] = vcf_snps(vcf)
    clair = {}
    for vcf in sorted((H / "clair3").glob("*/snps.vcf.gz")):
        clair.update(vcf_snps(vcf, vcf.parent.name))
    if clair:
        out["clair3"] = clair
    return out


def read_aln(path: Path) -> dict:
    seqs: dict = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def alignments() -> dict:
    """{name: (alignment, positions)}: columns are reference positions (all, or those listed)"""
    out = {"clair3": (H / "compare" / "clair3.aln", None)}
    for run in RUNS:
        d = H / "bacon" / run / "4_compared" / "minimap2"
        if (d / "minimap2.positions.tsv").is_file():
            lines = (d / "minimap2.positions.tsv").read_text().splitlines()
            out[f"bacon_minimap2_{run}"] = (d / "minimap2.variable.fasta", [int(x.split("\t")[1]) for x in lines])
    return {k: v for k, v in out.items() if v[0].is_file()}


def summary(err: list[int]) -> list:
    n = len(err) or 1
    return [sum(e == 0 for e in err), sum(abs(e) <= 1 for e in err), sum(abs(e) <= 5 for e in err),
            round(sum(map(abs, err)) / n, 2), round(sum(err) / n, 2), max(map(abs, err), default=0)]


def main() -> None:
    truth = matrix(H / "truth" / "distances.tsv")
    names_t = {a for a, _ in truth}
    zero = ZERO + [(a, b) for (a, b), v in truth.items() if v == 0 and a < b and (a, b) not in ZERO]
    (H / "errors").mkdir(exist_ok=True)
    head = ["exact", "within_1", "within_5", "mean_abs_error", "mean_signed_error", "max_abs_error"]
    with open(H / "scores_distances.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["method", "samples", "pairs", *head, "zero_pairs"])
        for name, path in methods().items():
            m = matrix(path)
            samples = sorted({a for a, _ in m} & names_t)
            pairs = [(a, b) for i, a in enumerate(samples) for b in samples[i + 1:]]
            err = [m[p] - truth[p] for p in pairs]
            with open(H / "errors" / f"{name}.tsv", "w") as out:
                out.write("a\tb\ttruth\tmethod\terror\n")
                out.writelines(f"{a}\t{b}\t{truth[a, b]}\t{m[a, b]}\t{e}\n" for (a, b), e in zip(pairs, err,
                                                                                               strict=True))
            zeros = "; ".join(f"{a}-{b}: {m.get((a, b), 'NA')}" for a, b in zero)
            w.writerow([name, len(samples), len(pairs), *summary(err), zeros])
    samples_t, table = truth_sites()
    with open(H / "scores_distances_judgeable.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["method", "samples", "pairs", *head])
        for name, (aln, positions) in alignments().items():
            seqs = read_aln(aln)
            names = sorted(set(seqs) & set(samples_t))
            err = []
            for i, a in enumerate(names):
                for b in names[i + 1:]:
                    masked = {p for p, v in table.items() if v[1][a] == "N" or v[1][b] == "N"}
                    d = sum(1 for k, (x, y) in enumerate(zip(seqs[a], seqs[b], strict=True))
                            if x != y and x in "ACGT" and y in "ACGT"
                            and (positions[k] if positions else k + 1) not in masked)
                    err.append(d - truth[a, b])
            w.writerow([name, len(names), len(err), *summary(err)])
    with open(H / "scores_snps.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["method", "samples", "true_snps", "found", "recall", "calls_at_callable", "true_calls",
                    "precision"])
        for name, calls in snp_calls().items():
            found = true_snps = called = true_calls = 0
            samples = [s for s in calls if s in samples_t]
            for s in samples_t:  # A sample without any call still has its true SNPs
                t = {p: v[1][s] for p, v in table.items() if v[1][s] not in ("N", v[0])}
                uncallable = {p for p, v in table.items() if v[1][s] == "N"}
                c = calls.get(s, {})
                true_snps += len(t)
                found += sum(c.get(p) == b for p, b in t.items())
                at = {p: b for p, b in c.items() if p not in uncallable}
                called += len(at)
                true_calls += sum(t.get(p) == b for p, b in at.items())
            w.writerow([name, len(samples), true_snps, found, round(found / (true_snps or 1), 3), called,
                        true_calls, round(true_calls / (called or 1), 3)])


if __name__ == "__main__":
    main()
