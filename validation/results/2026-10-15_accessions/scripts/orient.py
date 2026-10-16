#!/usr/bin/env python3
"""Put each plastome on the reference's frame before a multiple alignment, as done by hand in plastome studies:
of ptGAUL's two paths (the two orientations of the small single-copy region) and their reverse complements, the
one on which the most 2 kb windows (every 1 kb) of the reference map on the plus strand (minimap2 asm5); rotated
so that the base aligned to the reference's first base comes first.

Usage: orient.py REFERENCE.fasta OUT.fasta NAME PATH.fasta [PATH.fasta ...]"""
import subprocess
import sys
import tempfile

COMP = str.maketrans("ACGTN", "TGCAN")


def read_one(path: str) -> str:
    with open(path) as fh:
        return "".join(line.strip() for line in fh if not line.startswith(">")).upper()


def paf(target: str, query: str) -> list[list[str]]:
    with tempfile.NamedTemporaryFile("w", suffix=".fa") as t, tempfile.NamedTemporaryFile("w", suffix=".fa") as q:
        t.write(f">t\n{target}\n")
        q.write(f">q\n{query}\n")
        t.flush()
        q.flush()
        out = subprocess.run(["minimap2", "-x", "asm5", "-c", t.name, q.name], capture_output=True, text=True,
                             check=True).stdout
    return [line.split("\t") for line in out.splitlines()]


def paf_multi(target: str, queries: str) -> list[list[str]]:
    with tempfile.NamedTemporaryFile("w", suffix=".fa") as t, tempfile.NamedTemporaryFile("w", suffix=".fa") as q:
        t.write(f">t\n{target}\n")
        q.write(queries)
        t.flush()
        q.flush()
        out = subprocess.run(["minimap2", "-x", "asm5", "-c", t.name, q.name], capture_output=True, text=True,
                             check=True).stdout
    return [line.split("\t") for line in out.splitlines()]


def main() -> None:
    reference, out, name, paths = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
    ref = read_one(reference)
    # 2 kb windows of the reference every 1 kb, mapped to each path and its reverse complement (doubled: a
    # circle): the arrangement with the most windows on the plus strand is the reference's
    windows = "".join(f">w{i}\n{ref[i:i + 2000]}\n" for i in range(0, len(ref) - 2000, 1000))
    candidates = []
    for p in paths:
        for seq in (read_one(p), read_one(p).translate(COMP)[::-1]):
            hits = paf_multi(seq + seq, windows)
            plus = {f[0] for f in hits if f[4] == "+"}
            candidates.append((len(plus), seq))
    if not candidates:
        sys.exit(f"{name}: no alignment to the reference")
    seq = max(candidates)[1]
    # The reference's first 2 kb on the doubled sequence (a circle): the plus-strand hit gives the rotation
    hits = [f for f in paf(seq + seq, ref[:2000]) if f[4] == "+"]
    start = 0
    if hits:
        f = max(hits, key=lambda f: int(f[10]))
        start = (int(f[7]) - int(f[2])) % len(seq)
    seq = seq[start:] + seq[:start]
    with open(out, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")


if __name__ == "__main__":
    main()
