#!/usr/bin/env python3
"""Simulate validation datasets with a known truth for BACoN.

Scenarios
---------
plastid  A 70 kb circular chloroplast-like genome (LSC 40 kb, IRa 10 kb, SSC 10 kb, IRb = reverse complement of
         IRa) and nine samples on a known tree:
             s01, s02   identical to the reference (and to each other)
             s03        1 SNP;  s04 = s03 + 1 SNP          (very close samples)
             s05        clade X: 15 SNPs
             s06        X + 5 SNPs
             s07        X + 5 other SNPs + 2 SNPs inside the inverted repeat (both copies) + a 5 bp insertion
             s08        40 SNPs + a 1,500 bp deletion in the LSC (core versus pan genome)
             s09        s06 + 1 SNP, at low depth (10x)
rdna     A linear 22 kb molecule: 2 kb flank, a 3 kb unit repeated 6 times (like an rDNA array), 2 kb flank.
         Samples: t01 identical; t02 3 flank SNPs + 1 SNP in every copy of the unit; t03 = t02 + 2 SNPs;
         t04 5 flank SNPs + another unit SNP. Tests collapsed or circularized arrays.
linear   Three linear molecules (15, 9 and 6 kb; reads never wrap around) and four samples:
             l01 identical to the reference; l02 5 SNPs; l03 = l02 + 3 SNPs; l04 12 SNPs

Reads mimic R10.4.1 SUP reads: ~0.4% substitutions, ~0.2% random indels and homopolymer length errors that grow
with the run length; read lengths are log-normal (N50 ~8 kb); 15% of the reads come from an unrelated sequence
(baiting must remove them). Deterministic: the same seed gives the same files. Standard library only.

Writes, for each scenario: reference.fasta, reads/<sample>.fastq.gz, truth/<sample>.fasta (true genomes) and
truth.json (SNP events, expected distances, circularity, repeats).

Usage: python simulate.py OUTPUT_FOLDER [--scenario plastid|linear|all] [--depth 40] [--seed N]
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
from dataclasses import dataclass, field
from pathlib import Path

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def random_sequence(rng: random.Random, length: int, gc: float = 0.37) -> str:
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights=weights, k=length))


# ---------------------------------------------------------------------------------------------------------------
# Genomes
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class Molecule:
    name: str
    seq: str
    circular: bool


@dataclass
class Variant:
    """A change relative to the reference. `positions` are 0-based reference coordinates on `molecule`
    (two positions for a SNP inside an inverted repeat: the change is applied to both copies)."""

    id: str
    kind: str  # snp, ins, del
    molecule: int
    positions: list[int]
    alt: str = ""  # SNP: new base at positions[0]; ins: inserted bases
    length: int = 0  # del: deleted length
    in_repeat: bool = False
    inverted: bool = True  # Copies after the first are reverse-complemented (inverted repeat) or not (tandem)


@dataclass
class Scenario:
    name: str
    reference: list[Molecule]
    samples: dict[str, list[str]] = field(default_factory=dict)  # sample -> variant ids
    variants: dict[str, Variant] = field(default_factory=dict)
    depth: dict[str, float] = field(default_factory=dict)
    repeats: list[dict] = field(default_factory=list)


def _choose_positions(rng: random.Random, n: int, allowed: list[tuple[int, int]], used: set[int]) -> list[int]:
    out = []
    while len(out) < n:
        start, end = rng.choice(allowed)
        p = rng.randrange(start, end)
        if all(abs(p - u) > 50 for u in used):  # Keep variants apart (and away from each other's k-mers)
            used.add(p)
            out.append(p)
    return out


def _alt(rng: random.Random, base: str) -> str:
    return rng.choice([b for b in "ACGT" if b != base])


def plastid(rng: random.Random) -> Scenario:
    lsc, ir, ssc = 40_000, 10_000, 10_000
    lsc_seq, ir_seq, ssc_seq = random_sequence(rng, lsc), random_sequence(rng, ir), random_sequence(rng, ssc)
    genome = lsc_seq + ir_seq + ssc_seq + revcomp(ir_seq)
    irb_start = lsc + ir + ssc
    sc = Scenario("plastid", [Molecule("plastid", genome, True)])
    sc.repeats = [{"name": "IRa", "start": lsc, "end": lsc + ir}, {"name": "IRb", "start": irb_start,
                                                                     "end": irb_start + ir}]
    single_copy = [(300, lsc - 300), (lsc + ir + 300, irb_start - 300)]
    used: set[int] = set()

    def snps(prefix: str, n: int) -> list[str]:
        ids = []
        for i, p in enumerate(_choose_positions(rng, n, single_copy, used)):
            vid = f"{prefix}{i + 1}"
            sc.variants[vid] = Variant(vid, "snp", 0, [p], _alt(rng, genome[p]))
            ids.append(vid)
        return ids

    x = snps("X", 15)
    s03 = snps("a", 1)
    s04 = snps("b", 1)
    s06 = snps("c", 5)
    s07 = snps("d", 5)
    s08 = snps("e", 40)
    s09 = snps("f", 1)
    # Two SNPs in IRa, mirrored in IRb.
    ir_ids = []
    for i, off in enumerate(_choose_positions(rng, 2, [(lsc + 500, lsc + ir - 500)], used)):
        o = off - lsc
        mirror = irb_start + (ir - 1 - o)
        vid = f"IR{i + 1}"
        sc.variants[vid] = Variant(vid, "snp", 0, [off, mirror], _alt(rng, genome[off]), in_repeat=True)
        ir_ids.append(vid)
    ins_pos = _choose_positions(rng, 1, single_copy, used)[0]
    sc.variants["ins1"] = Variant("ins1", "ins", 0, [ins_pos], random_sequence(rng, 5))
    # A deletion in the LSC, clear of every other variant.
    while True:
        start = rng.randrange(2_000, lsc - 4_000)
        if all(not (start - 100 <= u <= start + 1_600) for u in used):
            break
    sc.variants["del1"] = Variant("del1", "del", 0, [start], length=1_500)
    sc.samples = {
        "s01": [], "s02": [], "s03": s03, "s04": s03 + s04, "s05": x, "s06": x + s06,
        "s07": x + s07 + ir_ids + ["ins1"], "s08": s08 + ["del1"], "s09": x + s06 + s09,
    }
    sc.depth = {"s09": 10.0}
    return sc


def linear(rng: random.Random) -> Scenario:
    sizes = (15_000, 9_000, 6_000)
    mols = [Molecule(f"linear{i + 1}", random_sequence(rng, n), False) for i, n in enumerate(sizes)]
    sc = Scenario("linear", mols)
    used: dict[int, set[int]] = {i: set() for i in range(len(mols))}

    def snps(prefix: str, n: int) -> list[str]:
        ids = []
        for i in range(n):
            m = rng.randrange(len(mols))
            p = _choose_positions(rng, 1, [(300, len(mols[m].seq) - 300)], used[m])[0]
            vid = f"{prefix}{i + 1}"
            sc.variants[vid] = Variant(vid, "snp", m, [p], _alt(rng, mols[m].seq[p]))
            ids.append(vid)
        return ids

    l02 = snps("g", 5)
    l03 = snps("h", 3)
    l04 = snps("k", 12)
    sc.samples = {"l01": [], "l02": l02, "l03": l02 + l03, "l04": l04}
    return sc


def rdna(rng: random.Random) -> Scenario:
    flank, unit, copies = 2_000, 3_000, 6
    unit_seq = random_sequence(rng, unit, gc=0.55)
    left, right = random_sequence(rng, flank), random_sequence(rng, flank)
    genome = left + unit_seq * copies + right
    sc = Scenario("rdna", [Molecule("rdna", genome, False)])
    sc.repeats = [{"name": "tandem_array", "start": flank, "end": flank + unit * copies, "unit": unit,
                   "copies": copies}]
    used: set[int] = set()
    flanks = [(300, flank - 300), (flank + unit * copies + 300, len(genome) - 300)]

    def snps(prefix: str, n: int) -> list[str]:
        ids = []
        for i, p in enumerate(_choose_positions(rng, n, flanks, used)):
            vid = f"{prefix}{i + 1}"
            sc.variants[vid] = Variant(vid, "snp", 0, [p], _alt(rng, genome[p]))
            ids.append(vid)
        return ids

    def unit_snp(vid: str) -> str:
        off = rng.randrange(200, unit - 200)
        positions = [flank + c * unit + off for c in range(copies)]
        sc.variants[vid] = Variant(vid, "snp", 0, positions, _alt(rng, unit_seq[off]), in_repeat=True,
                                   inverted=False)
        return vid

    t02 = snps("m", 3) + [unit_snp("U1")]
    t03 = snps("n", 2)
    t04 = snps("p", 5) + [unit_snp("U2")]
    sc.samples = {"t01": [], "t02": t02, "t03": t02 + t03, "t04": t04}
    return sc


def apply_variants(sc: Scenario, variant_ids: list[str]) -> list[Molecule]:
    """The sample's genome: the reference with its variants applied (right to left, so coordinates hold)."""
    seqs = [list(m.seq) for m in sc.reference]
    edits = []
    for vid in variant_ids:
        v = sc.variants[vid]
        if v.kind == "snp":
            for k, p in enumerate(v.positions):
                # The copy in an inverted repeat is reverse-complemented; tandem copies are identical.
                base = v.alt if k == 0 or not v.inverted else v.alt.translate(COMPLEMENT)
                edits.append((v.molecule, p, "snp", base, 0))
        else:
            edits.append((v.molecule, v.positions[0], v.kind, v.alt, v.length))
    for mol, pos, kind, alt, length in sorted(edits, key=lambda e: (e[0], e[1]), reverse=True):
        s = seqs[mol]
        if kind == "snp":
            s[pos] = alt
        elif kind == "ins":
            s[pos + 1: pos + 1] = list(alt)
        else:
            del s[pos: pos + length]
    return [Molecule(m.name, "".join(s), m.circular) for m, s in zip(sc.reference, seqs)]


# ---------------------------------------------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------------------------------------------

SUB, INDEL = 0.004, 0.002


def add_errors(seq: str, rng: random.Random) -> str:
    out = []
    i, n = 0, len(seq)
    while i < n:
        j = i
        while j < n and seq[j] == seq[i]:
            j += 1
        run = j - i
        base = seq[i]
        # Homopolymer length errors: more likely in long runs, mostly deletions (typical of R10.4.1 SUP).
        if run >= 3 and rng.random() < 0.012 * (run - 2):
            run += 1 if rng.random() < 0.35 else -1
        for _ in range(run):
            r = rng.random()
            if r < SUB:
                out.append(_alt(rng, base))
            elif r < SUB + INDEL / 2:
                continue
            elif r < SUB + INDEL:
                out.append(base + rng.choice("ACGT"))
            else:
                out.append(base)
        i = j
    return "".join(out)


def read_length(rng: random.Random) -> int:
    return max(800, min(40_000, int(rng.lognormvariate(8.7, 0.6))))


def quality(rng: random.Random, n: int) -> str:
    return "".join(chr(33 + max(3, min(50, int(rng.gauss(22, 6))))) for _ in range(n))


def simulate_reads(genome: list[Molecule], off_target: str, depth: float, rng: random.Random, name: str,
                   off_target_fraction: float = 0.15) -> list[str]:
    total = sum(len(m.seq) for m in genome)
    weights = [len(m.seq) for m in genome]
    records, bases, n = [], 0, 0
    while bases < depth * total:
        n += 1
        length = read_length(rng)
        if rng.random() < off_target_fraction:
            start = rng.randrange(len(off_target) - length)
            fragment = off_target[start: start + length]
        else:
            mol = rng.choices(genome, weights=weights)[0]
            if mol.circular:
                start = rng.randrange(len(mol.seq))
                fragment = (mol.seq + mol.seq)[start: start + min(length, len(mol.seq))]
            else:  # Linear molecule: fragments stop at the ends
                start = rng.randrange(-length + 200, len(mol.seq) - 200)
                fragment = mol.seq[max(0, start): max(0, start) + length if start >= 0 else start + length]
            if len(fragment) < 200:
                continue
            bases += len(fragment)
        if rng.random() < 0.5:
            fragment = revcomp(fragment)
        seq = add_errors(fragment, rng)
        records.append(f"@{name}_{n} simulated\n{seq}\n+\n{quality(rng, len(seq))}\n")
    return records


# ---------------------------------------------------------------------------------------------------------------
# Truth
# ---------------------------------------------------------------------------------------------------------------

def expected_distances(sc: Scenario) -> dict[str, dict[str, dict[str, int]]]:
    """Pairwise SNP counts under three conventions:
    events   each SNP counted once (a SNP in the inverted repeat is one event)
    sites    each changed position (a SNP in the inverted repeat counts twice)
    no_repeat SNPs outside repeats only (what reference mapping reports when repeats get MAPQ 0)
    Indels and the deletion are not SNPs and are never counted.
    """
    names = ["Reference", *sc.samples]
    snps = {n: {v for v in sc.samples.get(n, []) if sc.variants[v].kind == "snp"} for n in names}
    out: dict[str, dict[str, dict[str, int]]] = {"events": {}, "sites": {}, "no_repeat": {}}
    for a in names:
        for conv in out:
            out[conv][a] = {}
        for b in names:
            diff = snps[a] ^ snps[b]
            out["events"][a][b] = len(diff)
            out["sites"][a][b] = sum(len(sc.variants[v].positions) for v in diff)
            out["no_repeat"][a][b] = sum(1 for v in diff if not sc.variants[v].in_repeat)
    return out


def write_fasta(path: Path, molecules: list[Molecule]) -> None:
    with open(path, "w") as fh:
        for m in molecules:
            fh.write(f">{m.name} {'circular' if m.circular else 'linear'}\n")
            for i in range(0, len(m.seq), 80):
                fh.write(m.seq[i: i + 80] + "\n")


def build(sc: Scenario, out: Path, depth: float, rng: random.Random) -> None:
    (out / "reads").mkdir(parents=True, exist_ok=True)
    (out / "truth").mkdir(exist_ok=True)
    write_fasta(out / "reference.fasta", sc.reference)
    off_target = random_sequence(rng, 300_000, gc=0.45)
    for sample, variant_ids in sc.samples.items():
        genome = apply_variants(sc, variant_ids)
        write_fasta(out / "truth" / f"{sample}.fasta", genome)
        reads = simulate_reads(genome, off_target, sc.depth.get(sample, depth), rng, sample)
        # mtime=0: the same seed gives byte-identical files.
        with gzip.GzipFile(out / "reads" / f"{sample}.fastq.gz", "wb", compresslevel=6, mtime=0) as fh:
            fh.write("".join(reads).encode())
    truth = {
        "scenario": sc.name,
        "molecules": [{"name": m.name, "length": len(m.seq), "circular": m.circular} for m in sc.reference],
        "repeats": sc.repeats,
        "samples": {s: {"variants": v, "depth": sc.depth.get(s, depth)} for s, v in sc.samples.items()},
        "variants": {k: v.__dict__ for k, v in sc.variants.items()},
        "distances": expected_distances(sc),
    }
    (out / "truth.json").write_text(json.dumps(truth, indent=1) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("output", type=Path)
    parser.add_argument("--scenario", choices=["plastid", "linear", "rdna", "all"], default="all")
    parser.add_argument("--depth", type=float, default=40.0)
    parser.add_argument("--seed", type=int, default=20260924)
    args = parser.parse_args()
    for name, make in (("plastid", plastid), ("linear", linear), ("rdna", rdna)):
        if args.scenario in (name, "all"):
            rng = random.Random(f"{args.seed}-{name}")
            build(make(rng), args.output / name, args.depth, rng)
            print(f"{name}: written to {args.output / name}")


if __name__ == "__main__":
    main()
