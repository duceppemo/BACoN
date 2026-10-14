#!/usr/bin/env python3
"""Tables of the benchmark record from the scored runs (score_all.sh) and the /usr/bin/time files of the runs:
scores.tsv (per run and line, with the contigs and length of each assembly), distances.tsv (per run: SNPs to the
reference and the largest between two lines), times.tsv (wall time and peak memory) and
erratum_truth_covered.tsv (the 2026-10-13_mito record's truth_covered and the value rescored here).

Usage: BENCH=... summarize.py OUTPUT_FOLDER"""
import csv
import os
import sys
from pathlib import Path

B = Path(os.environ.get("BENCH", ".")) / "col0"
C = B / "compare"
MITO = Path(__file__).resolve().parents[2] / "2026-10-13_mito" / "col0_scores.tsv"
# (tool, organelle, folder of /usr/bin/time files relative to compare/, file name pattern, threads, what else ran
# on the machine, what one time covers)
A, S = "alone", "shared with the other pipelines' runs"
TIMES = [
    ("BACoN, templated", "cp", "bacon_timing", "cp_samtools.time", 16, A, "six lines, one at a time"),
    ("BACoN, Flye", "cp", "bacon_timing", "cp_flye.time", 16, A, "six lines, one at a time"),
    ("BACoN, templated, --filter-by-reference", "mt", "bacon_timing", "mt_samtools_fref.time", 16, A,
     "six lines, one at a time"),
    ("BACoN, Flye, --filter-by-reference", "mt", "bacon_timing", "mt_flye_fref.time", 16, A,
     "six lines, one at a time"),
    ("ptGAUL", "cp", "ptgaul", "MA32_*.time", 32, S, "one line"),
    ("TIPPo", "cp+mt", "tippo_ont", "MA32_*.time", 16, S + ", three lines at a time", "one line"),
    ("Oatk, k 501", "cp+mt", "oatk", "MA32_*_k501.time", 16, "three lines at a time", "one line"),
    ("CLAW (Flye and Unicycler)", "cp", "claw", "claw.time", 16, S, "six lines"),
    # MA32_119 alone on the machine (run_solo.sh)
    ("ptGAUL", "cp", "solo", "ptgaul.time", 16, A, "MA32_119"),
    ("TIPPo", "cp+mt", "solo", "tippo.time", 16, A, "MA32_119"),
    ("Oatk, k 501", "cp+mt", "solo", "oatk.time", 16, A, "MA32_119"),
    ("CLAW (Flye and Unicycler)", "cp", "solo", "claw.time", 16, A, "MA32_119"),
]


def tsv(path: Path) -> list[dict]:
    return list(csv.DictReader(path.read_text().splitlines(), delimiter="\t"))


def fasta_stats(path: Path) -> tuple[int, int]:
    contigs = length = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                contigs += 1
            else:
                length += len(line.strip())
    return contigs, length


def assemblies(run: str) -> Path:
    if run.startswith("bacon_"):
        return B / "bacon" / run.removeprefix("bacon_") / "3_assembled" / "all_assemblies"
    return C / "collect" / run / "3_assembled" / "all_assemblies"


def main() -> None:
    out = Path(sys.argv[1])
    with open(out / "scores.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "sample", "substitutions", "indels", "N_vs_truth", "truth_covered", "contigs",
                    "assembly_length"])
        for f in sorted((C / "scores").glob("*.tsv")):
            for row in tsv(f):
                contigs, length = fasta_stats(assemblies(f.stem) / f"{row['sample']}.fasta")
                w.writerow([f.stem, *row.values(), contigs, length])
    with open(out / "distances.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "lines", "snps_to_reference_min", "snps_to_reference_max", "max_snps_between_lines"])
        for f in sorted((C / "distances").glob("*.tsv")):
            rows = [line.split("\t") for line in f.read_text().splitlines()]
            names, matrix = rows[0][1:], [[int(x) for x in r[1:]] for r in rows[1:]]
            i, n = names.index("Reference"), len(names)
            to_ref = [matrix[i][j] for j in range(n) if j != i]
            between = [matrix[a][b] for a in range(n) for b in range(a + 1, n) if i not in (a, b)]
            w.writerow([f.stem, n - 1, min(to_ref), max(to_ref), max(between)])
    with open(out / "times.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["tool", "organelle", "run", "threads", "machine", "covers", "wall_s", "peak_rss_GB"])
        for tool, organelle, folder, pattern, threads, machine, covers in TIMES:
            for f in sorted((C / folder).glob(pattern)):
                if not f.stat().st_size:  # Still running
                    continue
                # The last line: "<seconds> s, <kB> kB" (a line above it when the command exited non-zero)
                wall, rss = f.read_text().strip().splitlines()[-1].replace(" s,", "").replace(" kB", "").split()
                # Oatk's MA32_99 at k 501 is the trial's run (run_oatk_trial.sh), during the other pipelines' runs
                shared = S if f.name == "MA32_99_k501.time" else machine
                w.writerow([tool, organelle, f.name.removesuffix(".time"), threads, shared, covers,
                            round(float(wall)), round(int(rss) / 1e6, 1)])

    with open(out / "erratum_truth_covered.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "sample", "truth_covered_one_to_one", "truth_covered_1kb"])
        for row in tsv(MITO):
            run = "bacon_" + row["run"].replace(" --filter-by-reference", "_fref")
            new = {r["sample"]: r for r in tsv(C / "scores" / f"{run}.tsv")}
            w.writerow([row["run"], row["sample"], row["truth_covered"], new[row["sample"]]["truth_covered"]])


if __name__ == "__main__":
    main()
