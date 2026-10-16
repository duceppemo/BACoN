#!/usr/bin/env python3
"""Prototype: BACoN's templated consensus rebuilt from each sample's aligned.bam with samtools consensus -C CUTOFF
(BACoN's options otherwise), then the minimap2 comparison (bacon.compare.run_minimap2, core). Writes
OUT/assemblies/<sample>.fasta, OUT/n_bases.tsv, OUT/mm2/snp_distances.tsv. Usage: consc.py RUN_DIR OUT CUTOFF"""
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from bacon.compare import run_minimap2, snp_distances, write_distances
from bacon.seqio import Record, acgtn, read_records, write_fasta

run, out, cutoff = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
asm = out / "assemblies"
asm.mkdir(parents=True, exist_ok=True)
bams = sorted(run.glob("3_assembled/*/aligned.bam"))


def one(bam: Path) -> tuple[str, int]:
    name = bam.parent.name
    raw = asm / f".{name}.raw.fasta"
    subprocess.run(["samtools", "consensus", "-a", "-X", "r10.4_sup", "--show-ins", "yes", "--show-del", "no",
                    "-d", "3", "-C", cutoff, "-f", "fasta", "-o", str(raw), str(bam)], check=True,
                   stderr=subprocess.DEVNULL)
    recs = [Record(f"{name}_{r.name}", acgtn(r.seq)) for r in read_records(raw) if r.seq]
    raw.unlink()
    write_fasta(asm / f"{name}.fasta", recs)
    return name, sum(r.seq.count("N") for r in recs)


with ThreadPoolExecutor(8) as ex:
    ns = list(ex.map(one, bams))
(out / "n_bases.tsv").write_text("".join(f"{n}\t{c}\n" for n, c in ns))
logs = out / "logs"
logs.mkdir(exist_ok=True)
paths = {n: asm / f"{n}.fasta" for n, _ in ns}
_, snps = run_minimap2(run / "reference.fasta", paths, out / "mm2", logs, threads=8, min_freq=1)
names, matrix = snp_distances(list(read_records(snps)))
write_distances(out / "mm2" / "snp_distances.tsv", names, matrix)
print(out, "N median", sorted(c for _, c in ns)[len(ns) // 2], "sites", len(next(iter(read_records(snps))).seq))
