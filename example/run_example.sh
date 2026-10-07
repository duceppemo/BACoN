#!/usr/bin/env bash
# Generate the simulated example (an annotated 30 kb circular reference, four samples with known SNPs and a sample
# metadata file), run BACoN on it with the default settings, and check the results against the truth: the SNP
# distances, the region, gene and effect of every SNP, and the report's annotation and metadata.
#
# Usage: bash example/run_example.sh [OUTPUT_FOLDER] [THREADS]
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
out=${1:-example_output}
threads=${2:-4}

python "$here/make_example.py" "$out/data"
bacon -r "$out/data/reference.fasta" -i "$out/data/reads" -o "$out/bacon" -t "$threads" -p 2 \
    --annotation "$out/data/reference.gb" --metadata "$out/data/metadata.tsv"

python - "$out" "$here" <<'PY'
import sys
from pathlib import Path

sys.path.insert(0, sys.argv[2])
from make_example import EXPECTED_DISTANCES

out = Path(sys.argv[1])
compared = out / "bacon" / "4_compared" / "ska"  # The default SNP method
failures = []

# 1. The pairwise SNP distances
distances = compared / "snp_distances.tsv"
if not distances.exists():
    sys.exit(f"FAILED: no {distances}")
lines = distances.read_text().splitlines()
names = lines[0].split("\t")[1:]
dist = {row.split("\t")[0]: dict(zip(names, map(int, row.split("\t")[1:]))) for row in lines[1:]}
wrong = [f"{a}-{b}: {dist[a][b]} (expected {d})" for (a, b), d in EXPECTED_DISTANCES.items() if dist[a][b] != d]
if wrong:
    failures.append("SNP distances: " + "; ".join(wrong))
else:
    print(f"OK: all {len(EXPECTED_DISTANCES)} pairwise SNP distances match the truth ({distances})")

# 2. The region, gene, context and effect of every SNP, as BACoN annotates the run's VCF, against the truth
from bacon.annotation import annotate_snps, load_annotation  # noqa: E402
from bacon.report import read_vcf  # noqa: E402
from bacon.seqio import read_records  # noqa: E402

before = len(failures)
rows = [line.split("\t") for line in (out / "data" / "planted_effects.tsv").read_text().splitlines()]
expected = {int(r[0]): r[4:10] for r in rows[1:]}  # position -> region, gene, context, codon, amino acid, effect
snps, _ = read_vcf(compared / "snps.vcf")
positions = {s.pos for s in snps}
if positions != set(expected):
    failures.append(f"VCF positions: {sorted(positions)} (expected {sorted(expected)})")
records = list(read_records(out / "bacon" / "reference.fasta"))
annotation = load_annotation(out / "data" / "reference.gb", [(r.name, len(r.seq)) for r in records])
info = annotate_snps(annotation, [(s.chrom, s.pos, s.ref, s.alt) for s in snps], {r.name: r.seq for r in records})
for s in snps:
    a = info.get((s.chrom, s.pos))
    if a is None or s.pos not in expected:
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
for text, what in [('<table class="snps', "the SNP table"), ("in the LSC", "the regions"),
                   ("coloured by <b>group</b>", "the colour column"), ("site 2", "the metadata columns")]:
    if text not in report:
        failures.append(f"report.html lacks {what} ({text!r})")
if failures:
    sys.exit("FAILED: " + "\n        ".join(failures))
print(f"OK: report.html has the annotation and the metadata ({out / 'bacon' / 'report.html'})")
PY
