#!/usr/bin/env python3
"""Score plastome assemblies of the 24 accessions against their Illumina truth (GetOrganelle 1.7.7.1: the two
paths of a circular plastome, the two orientations of the small single-copy region, of which the better counts;
or the scaffolds of an incomplete one): nucmer + show-snps (MUMmer 4), as the Col-0 benchmark's score.py.
Per accession: substitutions (both bases ACGT), indel events, positions where the assembly has N, the share of the
truth covered by alignments of at least 1 kb.

Usage: BENCH=... assembly_score.py FOLDER_OF_FASTAS  (needs MUMmer 4)."""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

G = Path(os.environ.get("BENCH", ".")) / "acc" / "truth" / "getorganelle"


def truths(sample: str) -> list[Path]:
    return sorted((G / sample).glob("*path_sequence.fasta"))


def score(asm: Path, ref: Path, organelle: str) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        pre = f"{tmp}/x"
        subprocess.run(["nucmer", "--maxmatch", "-p", pre, str(ref), str(asm)], check=True, capture_output=True)
        with open(f"{pre}.1delta", "w") as out:
            subprocess.run(["delta-filter", "-1", f"{pre}.delta"], check=True, stdout=out)
        snps = subprocess.run(["show-snps", "-C", "-H", "-T", "-r", f"{pre}.1delta"], check=True,
                              capture_output=True, text=True).stdout
        # Coverage from every alignment of at least 1 kb, not only the one-to-one ones: those drop one copy of
        # each repeat (a plastome's inverted repeat), and so understated the coverage of complete genomes
        with open(f"{pre}.ldelta", "w") as out:
            subprocess.run(["delta-filter", "-l", "1000", f"{pre}.delta"], check=True, stdout=out)
        coords = subprocess.run(["show-coords", "-H", "-T", "-r", f"{pre}.ldelta"], check=True,
                                capture_output=True, text=True).stdout
    sub = n = 0
    indels = set()
    for line in snps.splitlines():
        f = line.split("\t")
        r, q = f[1], f[2]
        if r == "." or q == ".":
            indels.add((f[-2], int(f[0]), r == "."))
        elif q.upper() == "N":
            n += 1
        else:
            sub += 1
    spans = {}
    for line in coords.splitlines():
        f = line.split("\t")
        spans.setdefault(f[-2], []).append((int(f[0]), int(f[1])))
    lengths = lengths_of(ref)
    covered = {}
    for name, iv in spans.items():
        e = cov = 0
        for a, b in sorted(iv):
            cov += max(0, b - max(a - 1, e))
            e = max(e, b)
        covered[name] = cov
    # Every sequence of the truth file counts (one circular path, or the scaffolds of an incomplete plastome)
    share = sum(covered.values()) / sum(lengths.values())
    # Indel events: an insertion in the assembly is listed at one truth position per base, a deletion at
    # consecutive truth positions; each run is one event
    events, last = 0, None
    for seq, pos, insertion in sorted(indels):
        if not (last and last[0] == seq and last[2] == insertion and pos - last[1] <= (0 if insertion else 1)):
            events += 1
        last = (seq, pos, insertion)
    return {"substitutions": sub, "indels": events, "N_vs_truth": n, "truth_covered": round(share, 3)}


def lengths_of(fa: Path) -> dict:
    out, name = {}, None
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = 0
            else:
                out[name] += len(line.strip())
    return out


def main() -> None:
    folder = Path(sys.argv[1])
    print("sample\tsubstitutions\tindels\tN_vs_truth\ttruth_covered\ttruth")
    for asm in sorted(folder.glob("*.fasta")):
        best = None
        for ref in truths(asm.stem):
            try:
                s = score(asm, ref, "cp")
            except subprocess.CalledProcessError:  # MUMmer found no alignment at all
                continue
            # The better orientation: the higher share covered, then the fewer differences
            key = (-s["truth_covered"], s["substitutions"] + s["indels"])
            if best is None or key < best[0]:
                best = (key, s, "circular" if len(truths(asm.stem)) > 1 else "scaffolds")
        if best is None:
            print(f"{asm.stem}\tNA\tNA\tNA\t0.0\t" + ("none" if not truths(asm.stem) else "no alignment"))
        else:
            print(asm.stem, *best[1].values(), best[2], sep="\t")


if __name__ == "__main__":
    main()
