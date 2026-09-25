#!/usr/bin/env python3
"""Score BACoN outputs against the truth written by simulate.py.

Usage: python score.py SCENARIO_FOLDER RUN_FOLDER [RUN_FOLDER ...] [--markdown OUT.md] [--json OUT.json]

For each run (a BACoN output folder, labelled by its folder name):
  assembly  per sample: contigs, circular contigs versus true circular molecules, total length versus truth,
            fraction of the true genome covered, duplicated bases, and consensus errors (mismatches, inserted and
            deleted bases) from minimap2 alignments of the contigs to the true genome
  SNPs      per comparison method found in 4_compared/: samples present, and for each pair of samples the
            reported SNP distance versus the truth under three conventions (events, sites, no_repeat; see
            simulate.py); a pair is exact when it equals the truth under the convention that fits best overall
Needs minimap2 on PATH.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

CS_OP = re.compile(r"(:\d+|\*[a-z]{2}|[+\-][a-z]+|~[a-z]{2}\d+[a-z]{2})")


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def read_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line.strip()]


def merge_intervals(intervals: list[tuple[int, int]]) -> int:
    total, end = 0, -1
    for s, e in sorted(intervals):
        if s > end:
            total += e - s
            end = e
        elif e > end:
            total += e - end
            end = e
    return total


def score_assembly(assembly: Path, truth: Path, repeats: list[dict]) -> dict:
    """Consensus errors and completeness of an assembly against the true genome.

    Inverted repeats make placement ambiguous: target positions in the second copy are folded onto the first
    before measuring how much of the genome is covered. Duplicated bases: aligned bases beyond the genome's
    unique length plus the one expected repeat copy (for example the overlapping ends of an uncircularized
    contig, or a region assembled twice).
    """
    contigs = read_fasta(assembly)
    true = read_fasta(truth)
    true_len = sum(map(len, true.values()))
    ir = [(r["start"], r["end"]) for r in repeats if r["name"].startswith("IR")]
    fold = None
    if len(ir) == 2:
        (a0, a1), (b0, b1) = ir
        fold = (a0, a1, b0, b1)
        unique_len = true_len - (b1 - b0)
    else:
        unique_len = true_len

    def folded(start: int, end: int) -> list[tuple[int, int]]:
        if not fold:
            return [(start, end)]
        a0, a1, b0, b1 = fold
        out = []
        if start < b0:
            out.append((start, min(end, b0)))
        if end > b0:  # Inside IRb: mirror onto IRa
            s, e = max(start, b0), min(end, b1)
            out.append((a0 + (b1 - e), a0 + (b1 - s)))
        if end > b1:
            out.append((b1, end))
        return out

    proc = subprocess.run(["minimap2", "-cx", "asm5", "--cs", "--secondary=no", "-t", "4", str(truth),
                           str(assembly)], capture_output=True, text=True, check=True)
    target: dict[str, list[tuple[int, int]]] = {}
    query: dict[str, list[tuple[int, int]]] = {}
    aligned_t = mism = ins = dels = 0
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if "tp:A:P" not in f[12:]:
            continue
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), "")
        target.setdefault(f[5], []).extend(folded(int(f[7]), int(f[8])))
        query.setdefault(f[0], []).append((int(f[2]), int(f[3])))
        aligned_t += int(f[8]) - int(f[7])
        for op in CS_OP.findall(cs):
            if op[0] == "*":
                mism += 1
            elif op[0] == "+":
                ins += len(op) - 1
            elif op[0] == "-":
                dels += len(op) - 1
    covered = sum(merge_intervals(v) for v in target.values())
    query_unique = sum(merge_intervals(v) for v in query.values())
    asm_len = sum(map(len, contigs.values()))
    errors = mism + ins + dels
    return {
        "contigs": len(contigs),
        "length": asm_len,
        "true_length": true_len,
        "covered_fraction": round(covered / unique_len, 4),
        "duplicated_bases": max(0, aligned_t - covered - (true_len - unique_len)),
        "unaligned_bases": asm_len - query_unique,
        "mismatches": mism, "inserted": ins, "deleted": dels,
        "errors_per_100kb": round(errors / max(query_unique, 1) * 1e5, 1),
    }


def read_distances(path: Path) -> dict[str, dict[str, int]]:
    rows = read_tsv(path)
    names = list(rows[0].keys())[1:]
    first = list(rows[0].keys())[0]
    return {r[first]: {n: int(r[n]) for n in names} for r in rows}


def score_snps(dist: dict[str, dict[str, int]], truth: dict) -> dict:
    expected = truth["distances"]
    names = list(expected["events"])
    present = [n for n in names if n in dist]
    pairs = [(a, b) for i, a in enumerate(present) for b in present[i + 1:]]
    by_conv = {}
    for conv, mat in expected.items():
        errs = [abs(dist[a][b] - mat[a][b]) for a, b in pairs]
        by_conv[conv] = {"exact_pairs": sum(e == 0 for e in errs), "max_abs_error": max(errs, default=0),
                         "total_abs_error": sum(errs)}
    best = min(by_conv, key=lambda c: (by_conv[c]["total_abs_error"], c))
    wrong = [f"{a}-{b}: {dist[a][b]} (truth {expected[best][a][b]})" for a, b in pairs
             if dist[a][b] != expected[best][a][b]]
    return {"samples_expected": len(names), "samples_present": len(present),
            "missing": [n for n in names if n not in dist], "pairs": len(pairs), "best_convention": best,
            **by_conv[best], "wrong_pairs": wrong, "by_convention": by_conv}


def score_run(scenario: Path, run: Path, truth: dict) -> dict:
    result: dict = {"run": run.name, "assembly": {}, "snps": {}}
    summary = {r["Sample"]: r for r in read_tsv(run / "summary.tsv")} if (run / "summary.tsv").exists() else {}
    circular_truth = sum(m["circular"] for m in truth["molecules"])
    for sample in truth["samples"]:
        asm = run / "3_assembled" / "all_assemblies" / f"{sample}.fasta"
        row = summary.get(sample, {})
        if not asm.is_file():
            result["assembly"][sample] = {"status": row.get("Status", "missing")}
            continue
        scored = score_assembly(asm, scenario / "truth" / f"{sample}.fasta", truth["repeats"])
        scored["circular_contigs"] = row.get("Circular_contigs", "NA")
        scored["true_circular"] = circular_truth
        result["assembly"][sample] = scored
    for dist in sorted((run / "4_compared").glob("*/snp_distances.tsv")):
        result["snps"][dist.parent.name] = score_snps(read_distances(dist), truth)
    timing = run.with_name(run.name + ".time")  # Written by run_validation.sh: first line = the full run
    if timing.exists():
        result["duration_s"] = float(timing.read_text().split()[0])
    return result


def markdown(results: list[dict], truth: dict) -> str:
    out = [f"## Scenario `{truth['scenario']}`", ""]
    out += ["### Assemblies", "",
            "| Run | Sample | Contigs | Circular (truth) | Length / truth | Covered | Duplicated bp | "
            "Mismatches | Ins bp | Del bp | Errors /100 kb |",
            "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in results:
        for sample, a in r["assembly"].items():
            if "contigs" not in a:
                out.append(f"| {r['run']} | {sample} | {a['status']} | | | | | | | | |")
                continue
            out.append(f"| {r['run']} | {sample} | {a['contigs']} | "
                       f"{a['circular_contigs']} ({a['true_circular']}) "
                       f"| {a['length']:,} / {a['true_length']:,} | {a['covered_fraction']:.4f} | "
                       f"{a['duplicated_bases']:,} | {a['mismatches']} | {a['inserted']} | {a['deleted']} | "
                       f"{a['errors_per_100kb']} |")
    out += ["", "### SNP distances", "",
            "| Run | Method | Samples | Best convention | Exact pairs | Max error | Wrong pairs |",
            "|---|---|---|---|---|---|---|"]
    for r in results:
        for method, s in r["snps"].items():
            wrong = "; ".join(s["wrong_pairs"][:6]) + (" …" if len(s["wrong_pairs"]) > 6 else "")
            missing = f" (missing: {', '.join(s['missing'])})" if s["missing"] else ""
            out.append(f"| {r['run']} | {method} | {s['samples_present']}/{s['samples_expected']}{missing} | "
                       f"{s['best_convention']} | {s['exact_pairs']}/{s['pairs']} | {s['max_abs_error']} | "
                       f"{wrong} |")
    out += ["", "| Run | Wall time of the full run (s) |", "|---|---|"]
    out += [f"| {r['run']} | {r.get('duration_s', 'NA')} |" for r in results]
    return "\n".join(out) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("scenario", type=Path)
    parser.add_argument("runs", type=Path, nargs="+")
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    truth = json.loads((args.scenario / "truth.json").read_text())
    results = [score_run(args.scenario, run, truth) for run in args.runs if run.is_dir()]
    md = markdown(results, truth)
    if args.markdown:
        args.markdown.write_text(md)
    if args.json:
        args.json.write_text(json.dumps(results, indent=1) + "\n")
    sys.stdout.write(md)


if __name__ == "__main__":
    main()
