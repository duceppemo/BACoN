#!/usr/bin/env python3
"""Masking test 1 (acc_masking_tests.tsv: w, wi and d variants). In each genome's projection on the reference
(bacon.compare.project, run before its N rule existed: rerun now, project applies that rule first), a base becomes
N when:

- w<W>: an N lies within W positions of it;
- w<W>i: an N, a deleted reference base or an insertion lies within W positions of it;
- d<W>k<K>: at least K N lie within W positions of it (the rule kept, with W 50 and K 3);

then core distances (columns with a nucleotide in every row, at least two different). Projections are cached as
<prefix><organelle>.proj.json. Usage: window.py OUT_PREFIX RUN_DIR cp|mt"""
import bisect
import json
import re
import subprocess
import sys
from pathlib import Path

from bacon.compare import project
from bacon.seqio import read_records

CS = re.compile(r"([=:*+\-~])([A-Za-z0-9]+)")


def projections(run: Path, cache: Path) -> dict:
    """{name: [projected sequence, offsets after which the genome has an insertion]}, the reference first"""
    if cache.is_file():
        return json.loads(cache.read_text())
    ref = list(read_records(run / "reference.fasta"))
    lengths = {r.name: len(r.seq) for r in ref}
    offsets, off = {}, 0
    for r in ref:
        offsets[r.name] = off
        off += len(r.seq)
    data = {"Reference": ["".join(r.seq.upper() for r in ref), []]}
    for p in sorted((run / "3_assembled" / "all_assemblies").glob("*.fasta")):
        cmd = ["minimap2", "-x", "asm20", "-c", "--cs=long", "-t", "16", str(run / "reference.fasta"), str(p)]
        paf = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        proj = project(paf, lengths)
        ins = []
        for line in paf.splitlines():
            f = line.split("\t")
            cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
            t = int(f[7]) + offsets[f[5]]
            for op, v in CS.findall(cs):
                if op == "+":
                    ins.append(t)
                t += len(v) if op in "=-" else 1 if op == "*" else 0
        data[p.stem] = ["".join("".join(proj[r.name]) for r in ref), ins]
    cache.write_text(json.dumps(data))
    return data


def near(seq: str, ins: list[int], w: int, indel: bool) -> str:
    """w and wi variants"""
    s = list(seq)
    bad = [i for i, c in enumerate(seq) if c == "N" or (indel and c == "-")] + (ins if indel else [])
    for i in bad:
        for j in range(max(0, i - w), min(len(s), i + w + 1)):
            if s[j] != "-":
                s[j] = "N"
    return "".join(s)


def dense(seq: str, w: int, k: int) -> str:
    """d variants"""
    ns = [i for i, c in enumerate(seq) if c == "N"]
    return "".join("N" if c in "ACGT" and bisect.bisect_right(ns, i + w) - bisect.bisect_left(ns, i - w) >= k
                   else c for i, c in enumerate(seq))


def main() -> None:
    prefix, run, organelle = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    data = projections(run, prefix.parent / f"{prefix.name}{organelle}.proj.json")
    names = list(data)
    variants = [(f"w{w}", lambda s, i, w=w: near(s, i, w, False)) for w in (0, 3, 5, 10, 15)]
    variants += [(f"w{w}i", lambda s, i, w=w: near(s, i, w, True)) for w in (3, 5, 10, 15)]
    variants += [(f"d{w}k{k}", lambda s, i, w=w, k=k: dense(s, w, k)) for w in (10, 20, 50) for k in (2, 3)]
    for tag, mask in variants:
        seqs = {n: mask(*data[n]) for n in names}
        cols = [i for i in range(len(seqs["Reference"])) if all(seqs[n][i] in "ACGT" for n in names)
                and len({seqs[n][i] for n in names}) > 1]
        rows = ["snp-dists\t" + "\t".join(names)]
        for a in names:
            rows.append(a + "\t" + "\t".join(str(sum(seqs[a][i] != seqs[b][i] for i in cols)) for b in names))
        (prefix.parent / f"{prefix.name}{tag}_{organelle}.tsv").write_text("\n".join(rows) + "\n")
        print(f"{prefix.name}{organelle}\t{tag}\tSNP sites {len(cols)}", flush=True)


if __name__ == "__main__":
    main()
