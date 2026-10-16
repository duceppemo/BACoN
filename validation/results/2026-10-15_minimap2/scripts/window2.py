#!/usr/bin/env python3
"""Masking test 2 (acc_masking_tests.tsv: e variants; and the potato and Col-0 checks): on top of project()
(with its N rule), a base becomes N when at least K indel events (an insertion,
or a run of deleted reference bases) lie within W positions of it. Core distances. Projection caches as window.py.
Usage: window2.py OUT_PREFIX RUN_DIR cp|mt"""
import bisect
import json
import re
import subprocess
import sys
from pathlib import Path

from bacon.compare import project
from bacon.seqio import read_records

prefix, run, org = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
ref = list(read_records(run / "reference.fasta"))
lengths = {r.name: len(r.seq) for r in ref}
cache = prefix.parent / f"{prefix.name}{org}.proj2.json"
if cache.is_file():
    data = json.loads(cache.read_text())
else:
    offsets, off = {}, 0
    for r in ref:
        offsets[r.name] = off
        off += len(r.seq)
    data = {"Reference": ["".join(r.seq.upper() for r in ref), []]}
    for p in sorted((run / "3_assembled" / "all_assemblies").glob("*.fasta")):
        cmd = ["minimap2", "-x", "asm20", "-c", "--cs=long", "-t", "16", str(run / "reference.fasta"), str(p)]
        paf = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        proj = project(paf, lengths)
        events = []
        for line in paf.splitlines():
            f = line.split("\t")
            cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
            t = int(f[7]) + offsets[f[5]]
            for op, v in re.findall(r"([=:*+\-~])([A-Za-z0-9]+)", cs):
                if op in "+-":
                    events.append(t)
                t += len(v) if op in "=-" else 1 if op == "*" else 0
        data[p.stem] = ["".join("".join(proj[r.name]) for r in ref), sorted(set(events))]
    cache.write_text(json.dumps(data))


def dense(seq: str, ev: list[int], w: int, k: int) -> str:
    if k == 0:
        return seq
    return "".join("N" if c in "ACGT" and bisect.bisect_right(ev, i + w) - bisect.bisect_left(ev, i - w) >= k
                   else c for i, c in enumerate(seq))


names = list(data)
for w, k in [(0, 0), (50, 3), (50, 4), (50, 5), (100, 4), (100, 6), (100, 8)]:
    seqs = {n: dense(data[n][0], data[n][1], w, k) for n in names}
    cols = [i for i in range(len(seqs["Reference"])) if all(seqs[n][i] in "ACGT" for n in names)
            and len({seqs[n][i] for n in names}) > 1]
    tag = f"e{w}k{k}"
    rows = ["snp-dists\t" + "\t".join(names)]
    for a in names:
        rows.append(a + "\t" + "\t".join(str(sum(seqs[a][i] != seqs[b][i] for i in cols)) for b in names))
    (prefix.parent / f"{prefix.name}{tag}_{org}.tsv").write_text("\n".join(rows) + "\n")
    masked = sum(s.count("N") for n, s in seqs.items()) - sum(data[n][0].count("N") for n in names)
    print(f"{prefix.name}{org}\t{tag}\tSNP sites {len(cols)}\tmasked {masked}", flush=True)
