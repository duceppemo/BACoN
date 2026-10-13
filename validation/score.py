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
heteroplasmy  per sample (a mixture of a major and a minor haplotype): what the assembly has at each site where
            the minor haplotype differs: the major allele, the minor allele, N, another sequence, both alleles in
            different contigs (split), or nothing (missing); and the reported SNP distance to the reference
numt        per sample: reads by origin (organelle, NUPT, nuclear; from the read names) among the simulated, the
            baited (1_extracted/) and the filtered (2_filtered/) reads; consensus errors inside the plastome
            segments copied into the nucleus and elsewhere; what the assembly has at the substitutions that
            distinguish the NUPTs from the plastome (the NUPT allele must never appear)
Needs minimap2 on PATH.
"""

from __future__ import annotations

import argparse
import gzip
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


def primary_alignments(assembly: Path, truth: Path) -> list[dict]:
    """Primary minimap2 alignments of the contigs to the true genome, with, for each target position they cover,
    the contig's bases there (in target orientation; '' where deleted) and the bases inserted after it."""
    true = read_fasta(truth)
    proc = subprocess.run(["minimap2", "-cx", "asm5", "--cs", "--secondary=no", "-t", "4", str(truth),
                           str(assembly)], capture_output=True, text=True, check=True)
    out = []
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if "tp:A:P" not in f[12:]:
            continue
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), "")
        start, end, seq = int(f[7]), int(f[8]), true[f[5]]
        bases, inserted, errors = [""] * (end - start), [""] * (end - start), []
        t = start
        for op in CS_OP.findall(cs):
            if op[0] == ":":
                n = int(op[1:])
                bases[t - start: t - start + n] = seq[t: t + n]
                t += n
            elif op[0] == "*":
                bases[t - start] = op[2].upper()
                errors.append((t, 1))
                t += 1
            elif op[0] == "+":
                if t > start:
                    inserted[t - 1 - start] += op[1:].upper()
                errors.append((t, len(op) - 1))
            elif op[0] == "-":
                errors.append((t, len(op) - 1))
                t += len(op) - 1
        out.append({"target": f[5], "start": start, "end": end, "bases": bases, "inserted": inserted,
                    "errors": errors})
    return out


def site_call(alignments: list[dict], target: str, windows: list[tuple[int, int, str, str]]) -> str:
    """What the assembly has over the target windows [a, b) (one, or a SNP and its mirror in the other copy of
    an inverted repeat): the major or the minor allele (each given as the whole window), N, other, split (the
    contigs, or the two repeat copies, disagree) or missing (no alignment covers a window)."""
    calls = set()
    for a, b, major, minor in windows:
        for aln in alignments:
            if aln["target"] != target or not aln["start"] <= a < b <= aln["end"]:
                continue
            o = aln["start"]
            q = "".join(aln["bases"][i - o] + (aln["inserted"][i - o] if i < b - 1 else "") for i in range(a, b))
            calls.add("major" if q == major else "minor" if q == minor else "N" if "N" in q else "other")
    return calls.pop() if len(calls) == 1 else "split:" + "+".join(sorted(calls)) if calls else "missing"


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def windows(site: dict, seq: str, repeats: list[dict], margin: int = 5) -> list[tuple[int, int, str, str]]:
    """The window [a, b) around a minor site and both alleles over it. Indels get a margin on each side, so that
    an aligner may place them anywhere in a repeat of their bases. A SNP in an inverted repeat also gets its
    mirror: a contig holds both copies, but either may be aligned to either copy of the reference."""
    p = site["position"]
    if site["kind"] == "snp":
        out = [(p, p + 1, site["major"], site["minor"])]
        ir = [(r["start"], r["end"]) for r in repeats if r["name"].startswith("IR")]
        if site.get("in_repeat") and len(ir) == 2:
            (a0, a1), (b0, b1) = ir
            m = b1 - 1 - (p - a0) if a0 <= p < a1 else a0 + (b1 - 1 - p)
            out.append((m, m + 1, site["major"].translate(COMPLEMENT), site["minor"].translate(COMPLEMENT)))
        return out
    if site["kind"] == "ins":
        a, b = p - margin, p + 1 + margin
        return [(a, b, seq[a:b], seq[a: p + 1] + site["minor"] + seq[p + 1: b])]
    a, b = p - margin, p + len(site["major"]) + margin
    return [(a, b, seq[a:b], seq[a:p] + seq[p + len(site["major"]): b])]


CALLS = ("major", "minor", "N", "other", "split", "missing")


def score_heteroplasmy(scenario: Path, run: Path, truth: dict, snps: dict) -> dict:
    out = {}
    dist = {m: read_distances(p) for m in snps for p in [run / "4_compared" / m / "snp_distances.tsv"]}
    for sample, t in truth["samples"].items():
        asm = run / "3_assembled" / "all_assemblies" / f"{sample}.fasta"
        row: dict = {"minor_fraction": t["minor_fraction"],
                     "distance_to_reference": {m: d.get("Reference", {}).get(sample) for m, d in dist.items()},
                     "true_distance": truth["distances"]["events"]["Reference"][sample]}
        if asm.is_file():
            target, seq = next(iter(read_fasta(scenario / "truth" / f"{sample}.fasta").items()))
            alns = primary_alignments(asm, scenario / "truth" / f"{sample}.fasta")
            sites = {}
            for site in t["minor_sites"]:
                sites[site["id"]] = site_call(alns, target, windows(site, seq, truth["repeats"]))
            row["sites"] = sites
            snp_calls = [sites[site["id"]].split(":")[0] for site in t["minor_sites"] if site["kind"] == "snp"]
            row["snp_calls"] = {c: snp_calls.count(c) for c in CALLS}
            row["indel_calls"] = {k: sites[k] for k in (x["id"] for x in t["minor_sites"] if x["kind"] != "snp")}
            row["n_bases"] = sum(s.count("N") for s in read_fasta(asm).values())
        out[sample] = row
    return out


def origins(path: Path) -> dict[str, list[int]]:
    """Reads and bases by origin, the last token of the read names written by simulate.py."""
    counts: dict[str, list[int]] = {}
    if not path.is_file():
        return counts
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                origin = line[1:].split()[0].rsplit("_", 1)[-1]
            elif i % 4 == 1:
                c = counts.setdefault(origin, [0, 0])
                c[0] += 1
                c[1] += len(line.strip())
    return counts


def score_numt(scenario: Path, run: Path, truth: dict, snps: dict, cache: dict) -> dict:
    out = {}
    regions = truth["copied_regions"]
    dist = {m: read_distances(p) for m in snps for p in [run / "4_compared" / m / "snp_distances.tsv"]}
    for sample, t in truth["samples"].items():
        if sample not in cache:
            cache[sample] = origins(scenario / "reads" / f"{sample}.fastq.gz")
        row: dict = {"nuclear_depth": t["nuclear_depth"], "nuclear_genome": t["nuclear_genome"],
                     "reads": {"simulated": cache[sample],
                               "baited": origins(run / "1_extracted" / f"{sample}.fastq.gz"),
                               "filtered": origins(run / "2_filtered" / f"{sample}.fastq.gz")},
                     "distance_to_reference": {m: d.get("Reference", {}).get(sample) for m, d in dist.items()},
                     "true_distance": truth["distances"]["events"]["Reference"][sample]}
        asm = run / "3_assembled" / "all_assemblies" / f"{sample}.fasta"
        if asm.is_file():
            target, seq = next(iter(read_fasta(scenario / "truth" / f"{sample}.fasta").items()))
            alns = primary_alignments(asm, scenario / "truth" / f"{sample}.fasta")
            inside = [0, 0]  # Error bp, aligned bp: inside the copied segments
            elsewhere = [0, 0]
            for aln in alns:
                within = sum(max(0, min(aln["end"], e) - max(aln["start"], s)) for s, e in regions)
                inside[1] += within
                elsewhere[1] += aln["end"] - aln["start"] - within
                for pos, bp in aln["errors"]:
                    (inside if any(s <= pos < e for s, e in regions) else elsewhere)[0] += bp
            row["errors_copied"], row["errors_elsewhere"] = inside, elsewhere
            calls: dict[str, dict[str, int]] = {}  # Insert -> call -> count, at the NUPT substitutions
            nupt_sites = []
            for ins in truth["nuclear"][t["nuclear_genome"]]["inserts"]:
                c = calls.setdefault(str(ins["id"]), dict.fromkeys(CALLS, 0))
                for pos, _, nupt in ins["substitutions"]:
                    if seq[pos] == nupt:  # The sample's own SNP gave the plastome the NUPT's base
                        continue
                    call = site_call(alns, target, [(pos, pos + 1, seq[pos], nupt)])
                    c[call.split(":")[0]] += 1
                    if call == "minor" or "minor" in call.split(":")[-1].split("+"):  # In any contig
                        nupt_sites.append(pos)
            row["nupt_sites"] = calls
            row["nupt_alleles_in_assembly"] = nupt_sites
            # The sample's own SNPs inside the copied segments, where the NUPT reads carry the reference base
            # (minor: the reference base was called).
            reference = next(iter(read_fasta(scenario / "reference.fasta").values()))
            own = {}
            for vid in t["variants"]:
                pos = truth["variants"][vid]["positions"][0]
                if any(s <= pos < e for s, e in regions):
                    own[vid] = site_call(alns, target, [(pos, pos + 1, seq[pos], reference[pos])])
            row["own_snps_copied"] = own
        out[sample] = row
    return out


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


def score_run(scenario: Path, run: Path, truth: dict, cache: dict | None = None) -> dict:
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
    if truth["scenario"] == "heteroplasmy":
        result["heteroplasmy"] = score_heteroplasmy(scenario, run, truth, result["snps"])
    elif truth["scenario"] == "numt":
        result["numt"] = score_numt(scenario, run, truth, result["snps"], {} if cache is None else cache)
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
    if truth["scenario"] == "heteroplasmy":
        out += markdown_heteroplasmy(results, truth)
    elif truth["scenario"] == "numt":
        out += markdown_numt(results, truth)
    out += ["", "| Run | Wall time of the full run (s) |", "|---|---|"]
    out += [f"| {r['run']} | {r.get('duration_s', 'NA')} |" for r in results]
    return "\n".join(out) + "\n"


def _distance(row: dict) -> str:
    reported = "; ".join(f"{m} {d}" for m, d in row["distance_to_reference"].items() if d is not None)
    return f"{reported or 'NA'} ({row['true_distance']})"


def markdown_heteroplasmy(results: list[dict], truth: dict) -> list[str]:
    n_snp = sum(s["kind"] == "snp" for s in next(iter(truth["samples"].values()))["minor_sites"])
    indels = [s["id"] for s in next(iter(truth["samples"].values()))["minor_sites"] if s["kind"] != "snp"]
    out = ["", "### Minor-haplotype sites", "",
           f"At the {n_snp} SNP sites where the minor haplotype differs (inverted-repeat SNPs count twice), what "
           "the assembly has: the major allele (correct), the minor allele, N, another base, both alleles in "
           "different contigs or repeat copies (split) or nothing (missing); then the "
           f"{len(indels)} indels ({', '.join(indels)}): M major, m minor, N, o other, s split, - missing. "
           "Distance: SNPs to the reference reported by each method (truth).", "",
           "| Run | Sample | Minor fraction | Major | Minor | N | Other | Split | Missing | Indels | N bases | "
           "Distance to reference (truth) |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    short = {"major": "M", "minor": "m", "N": "N", "other": "o", "split": "s", "missing": "-"}
    for r in results:
        rows = sorted(r["heteroplasmy"].items(), key=lambda kv: (kv[1]["minor_fraction"], kv[0]))
        for sample, h in rows:
            if "snp_calls" not in h:
                out.append(f"| {r['run']} | {sample} | {h['minor_fraction']:.2f} | no assembly | | | | | | | | |")
                continue
            c = h["snp_calls"]
            pct = [f"{c[k]} ({100 * c[k] / n_snp:.0f}%)" for k in ("major", "minor", "N")]
            out.append(f"| {r['run']} | {sample} | {h['minor_fraction']:.2f} | {' | '.join(pct)} | {c['other']} | "
                       f"{c['split']} | {c['missing']} | "
                       f"{'/'.join(short[v.split(':')[0]] for v in h['indel_calls'].values())} | {h['n_bases']} | "
                       f"{_distance(h)} |")
    return out


def _reads(counts: dict[str, list[int]], kind: str) -> int:
    if kind == "nupt":
        return sum(v[0] for k, v in counts.items() if k.startswith("nupt"))
    return counts.get(kind, [0, 0])[0]


def markdown_numt(results: list[dict], truth: dict) -> list[str]:
    inserts = {g: v["inserts"] for g, v in truth["nuclear"].items()}
    out = ["", "### Reads by origin", "",
           "Reads simulated / baited (1_extracted) / kept by the filter (2_filtered). NUPT reads: nuclear reads "
           "covering at least 100 bp of an insert; nuclear: the other nuclear reads. Share: of the filtered "
           "reads (bases).", "",
           "| Run | Sample | Nuclear depth (genome) | Organelle reads | NUPT reads | Nuclear reads | "
           "NUPT share of filtered reads (bases) |",
           "|---|---|---|---|---|---|---|"]
    for r in results:
        for sample, n in r["numt"].items():
            rd = n["reads"]
            cells = ["/".join(str(_reads(rd[k], kind)) for k in ("simulated", "baited", "filtered"))
                     for kind in ("organelle", "nupt", "nuclear")]
            filt = rd["filtered"]
            total = [sum(v[i] for v in filt.values()) for i in (0, 1)]
            nupt = [sum(v[i] for k, v in filt.items() if k.startswith("nupt")) for i in (0, 1)]
            share = (f"{100 * nupt[0] / total[0]:.1f}% ({100 * nupt[1] / total[1]:.1f}%)" if total[0] else "NA")
            out.append(f"| {r['run']} | {sample} | {n['nuclear_depth']:g}x ({n['nuclear_genome']}) | "
                       f"{' | '.join(cells)} | {share} |")
    out += ["", "### NUPT reads baited, by insert", "",
            "Baited / simulated reads of each insert (length, identity to the plastome in genome A; in genome B "
            "every insert is 99.9% identical).", "",
            "| Run | Sample | " + " | ".join(f"#{i['id']} {i['length'] / 1000:g} kb {100 * i['identity']:g}%"
                                            for i in inserts["A"]) + " |",
            "|---|---|" + "---|" * len(inserts["A"])]
    for r in results:
        for sample, n in r["numt"].items():
            rd = n["reads"]
            cells = [f"{rd['baited'].get(f'nupt{i}', [0])[0]}/{rd['simulated'].get(f'nupt{i}', [0])[0]}"
                     for i in range(1, len(inserts["A"]) + 1)]
            out.append(f"| {r['run']} | {sample} | {' | '.join(cells)} |")
    out += ["", "### Consensus", "",
            "Consensus errors (mismatched, inserted and deleted bp) inside the plastome segments copied into the "
            "nucleus and elsewhere, per aligned length; at the substitutions that distinguish the NUPTs from "
            "the plastome, what the assembly has (plastome base / NUPT base / N / other / split: different bases "
            "in different contigs / missing) and the number of sites where the NUPT base is in any contig; the "
            "sample's own SNPs inside the copied segments (where the NUPT reads carry the reference base): "
            "called / reference base / N / split or other; distance: SNPs to the reference reported by each "
            "method (truth).", "",
            "| Run | Sample | Errors in copied segments (aligned bp) | Errors elsewhere (aligned bp) | "
            "NUPT sites: plastome / NUPT / N / other / split / missing | NUPT base in a contig | "
            "Own SNPs in copied segments | Distance to reference (truth) |",
            "|---|---|---|---|---|---|---|---|"]
    for r in results:
        for sample, n in r["numt"].items():
            if "errors_copied" not in n:
                out.append(f"| {r['run']} | {sample} | no assembly | | | | | |")
                continue
            tot = {k: sum(c[k] for c in n["nupt_sites"].values()) for k in CALLS}
            own = [v.split(":")[0] for v in n["own_snps_copied"].values()]
            own_cells = (f"{own.count('major')} / {own.count('minor')} / {own.count('N')} / "
                         f"{len(own) - own.count('major') - own.count('minor') - own.count('N')}")
            out.append(f"| {r['run']} | {sample} | {n['errors_copied'][0]} ({n['errors_copied'][1]:,}) | "
                       f"{n['errors_elsewhere'][0]} ({n['errors_elsewhere'][1]:,}) | "
                       f"{tot['major']} / {tot['minor']} / {tot['N']} / {tot['other']} / {tot['split']} / "
                       f"{tot['missing']} | {len(n['nupt_alleles_in_assembly'])} | {own_cells} | "
                       f"{_distance(n)} |")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("scenario", type=Path)
    parser.add_argument("runs", type=Path, nargs="+")
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    truth = json.loads((args.scenario / "truth.json").read_text())
    cache: dict = {}  # Origins of the simulated reads, counted once for all runs
    results = [score_run(args.scenario, run, truth, cache) for run in args.runs if run.is_dir()]
    md = markdown(results, truth)
    if args.markdown:
        args.markdown.write_text(md)
    if args.json:
        args.json.write_text(json.dumps(results, indent=1) + "\n")
    sys.stdout.write(md)


if __name__ == "__main__":
    main()
