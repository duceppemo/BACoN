#!/usr/bin/env python3
"""Score BACoN's Col-0 assemblies against the HiFi truth (TIPPo): nucmer + show-snps (MUMmer 4).
Per sample: substitutions (both bases ACGT), indel events, positions where the assembly has N, and the share of
the truth covered by alignments (show-coords).

Usage: BENCH=... col0_score.py RUN_FOLDER|GENOME.fasta cp|mt  (needs MUMmer 4)."""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

B = Path(os.environ.get("BENCH", ".")) / "col0"  # BENCH: the folder with col0/truth_tippo/


def truth(sample: str, organelle: str) -> Path | None:
    d = next((B / "truth_tippo" / sample).glob("*.organelle"), None)
    if d is None:
        return None
    if organelle == "cp":
        # The final plastome (both orientations of the single-copy region), or the last round's when TIPPo's
        # last round already gave one circular contig
        hits = sorted(d.glob("*.edge_*.organelle.chloroplast.fasta"), key=lambda p: p.stat().st_mtime)
        hits = hits or sorted(d.glob("*.round*.organelle.chloroplast.fasta"), key=lambda p: p.stat().st_mtime)
        return hits[-1] if hits else None
    asm = next(d.glob("*mitochondrial.fasta.filter.fasta.flye/assembly.fasta"), None)
    if asm is None:  # Rerun by hand where TIPPo's mitochondrial Flye failed
        asm = next((B / "truth_tippo" / sample).glob("mt_flye/assembly.fasta"), None)
    return asm


def score(asm: Path, ref: Path, organelle: str) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        pre = f"{tmp}/x"
        subprocess.run(["nucmer", "--maxmatch", "-p", pre, str(ref), str(asm)], check=True, capture_output=True)
        with open(f"{pre}.1delta", "w") as out:
            subprocess.run(["delta-filter", "-1", f"{pre}.delta"], check=True, stdout=out)
        snps = subprocess.run(["show-snps", "-C", "-H", "-T", "-r", f"{pre}.1delta"], check=True,
                              capture_output=True, text=True).stdout
        coords = subprocess.run(["show-coords", "-H", "-T", "-r", f"{pre}.1delta"], check=True,
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
    # Plastome: TIPPo gives both orientations of the single-copy region, one per sequence: the best one counts.
    # Mitogenome: contigs; all count.
    share = (max(covered.get(k, 0) / v for k, v in lengths.items()) if organelle == "cp"
             else sum(covered.values()) / sum(lengths.values()))
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
    run, organelle = Path(sys.argv[1]), sys.argv[2]
    if run.suffix in (".fasta", ".fa"):  # A single genome (the reference) against every truth
        for d in sorted((B / "truth_tippo").glob("MA32_*")):
            if (ref := truth(d.name, organelle)) is not None:
                print(d.name, *score(run, ref, organelle).values(), sep="\t")
        return
    print("sample\tsubstitutions\tindels\tN_vs_truth\ttruth_covered")
    for asm in sorted((run / "3_assembled" / "all_assemblies").glob("*.fasta")):
        ref = truth(asm.stem, organelle)
        if ref is None:
            print(f"{asm.stem}\tno truth")
            continue
        print(asm.stem, *score(asm, ref, organelle).values(), sep="\t")


if __name__ == "__main__":
    main()
