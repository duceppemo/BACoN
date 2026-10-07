"""Check a BACoN run of the bundled example against the truth written by make_example.py: the pairwise SNP
distances, the region, gene, context and effect of every SNP of the run's VCF, and the report's annotation and
metadata. Used by run_example.sh; `check(output_folder)` returns the failures and prints one OK line per check.

Usage: python example/check_example.py OUTPUT_FOLDER [SNP_METHOD]
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_example import EXPECTED_DISTANCES  # noqa: E402

from bacon.annotation import annotate_snps, load_annotation  # noqa: E402
from bacon.report import read_vcf  # noqa: E402
from bacon.seqio import read_records  # noqa: E402

SEQUENCE = "organelle"  # The reference's one sequence, as make_example names it

REPORT_TEXTS = [('<table class="snps', "the SNP table"), ("in the LSC", "the regions"),
                ("coloured by <b>group</b>", "the colour column"), ("site 2", "the metadata columns")]


def check(out: Path, method: str = "ska") -> list[str]:
    """The failures of the three checks (empty when the run agrees with the truth); an OK line is printed for each
    check that passes."""
    out = Path(out)
    compared = out / "bacon" / "4_compared" / method
    failures: list[str] = []

    # 1. The pairwise SNP distances
    distances = compared / "snp_distances.tsv"
    if not distances.exists():
        return [f"no {distances}"]
    lines = distances.read_text().splitlines()
    names = lines[0].split("\t")[1:]
    dist = {row.split("\t")[0]: dict(zip(names, map(int, row.split("\t")[1:]))) for row in lines[1:]}
    wrong = [f"{a}-{b}: {dist[a][b]} (expected {d})" for (a, b), d in EXPECTED_DISTANCES.items()
             if dist.get(a, {}).get(b) != d]
    if wrong:
        failures.append("SNP distances: " + "; ".join(wrong))
    else:
        print(f"OK: all {len(EXPECTED_DISTANCES)} pairwise SNP distances match the truth ({distances})")

    # 2. The region, gene, context and effect of every SNP, as BACoN annotates the run's VCF, against the truth
    before = len(failures)
    rows = [line.split("\t") for line in (out / "data" / "planted_effects.tsv").read_text().splitlines()]
    expected = {int(r[0]): r[4:10] for r in rows[1:]}  # position -> region, gene, context, codon, aa, effect
    snps, _ = read_vcf(compared / "snps.vcf")
    sites = {(s.chrom, s.pos) for s in snps}
    if sites != {(SEQUENCE, p) for p in expected}:
        failures.append(f"VCF sites: {sorted(sites)} (expected {sorted((SEQUENCE, p) for p in expected)})")
    records = list(read_records(out / "bacon" / "reference.fasta"))
    annotation = load_annotation(out / "data" / "reference.gb", [(r.name, len(r.seq)) for r in records])
    info = annotate_snps(annotation, [(s.chrom, s.pos, s.ref, s.alt) for s in snps],
                         {r.name: r.seq for r in records})
    for s in snps:
        if s.chrom != SEQUENCE or s.pos not in expected:
            continue  # Reported above
        a = info.get((s.chrom, s.pos))
        if a is None:
            failures.append(f"SNP {s.pos}: not annotated (expected {expected[s.pos]})")
            continue
        found = [a.region, ", ".join(dict.fromkeys(g.name for g in a.genes)), a.context,
                 "; ".join(e.codons for e in a.effects), "; ".join(e.change for e in a.effects),
                 "; ".join(e.kind for e in a.effects)]
        if found != expected[s.pos]:
            failures.append(f"SNP {s.pos}: {found} (expected {expected[s.pos]})")
    if len(failures) == before:
        effects = [r[9] for r in rows[1:] if r[9]]
        print(f"OK: the {len(expected)} SNPs are annotated as planted ({len(effects)} coding changes: "
              + ", ".join(f"{effects.count(k)} {k}" for k in dict.fromkeys(effects)) + ")")

    # 3. The report shows the annotation (regions, SNP table) and the metadata (the colour column)
    report = (out / "bacon" / "report.html").read_text(encoding="utf-8")
    for text, what in REPORT_TEXTS:
        if text not in report:
            failures.append(f"report.html lacks {what} ({text!r})")
    if not failures:
        print(f"OK: report.html has the annotation and the metadata ({out / 'bacon' / 'report.html'})")
    return failures


def main(argv: list[str]) -> None:
    if len(argv) not in (2, 3):
        sys.exit("Usage: python check_example.py OUTPUT_FOLDER [SNP_METHOD]")
    failures = check(Path(argv[1]), argv[2] if len(argv) == 3 else "ska")
    if failures:
        sys.exit("FAILED: " + "\n        ".join(failures))


if __name__ == "__main__":
    main(sys.argv)
