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
heteroplasmy  A 70 kb plastome (as in plastid) and samples that are mixtures of two haplotypes: a major one and a
         minor one = major + 25 SNPs (15 LSC, 5 SSC, 5 in the inverted repeat, both copies) + 3 indels (1-5 bp).
         Minor fractions 0 (h000), 0.02 (h002), 0.05, 0.10, 0.20, 0.30, 0.40 and 0.50 (h050, a tie); h020s has a
         minor fraction of 0.20 and 12 SNPs of its own in both haplotypes, so that distances are not all zero.
         Reads are drawn from each haplotype in proportion, at 80x. The truth is the major haplotype.
numt     A 120 kb plastome (LSC 80 kb, IR 12 kb, SSC 16 kb) at 50x and, instead of the unrelated off-target reads,
         reads from a 765 kb nuclear sequence holding eight copies of plastome segments (NUPTs) of 0.5, 2, 5 and
         15 kb at 99.9, 99, 97 and 90% identity (substitutions and 1-3 bp indels), between 80 kb of random
         sequence. Nuclear depth: n01 0.5x, n02 2x, n03 10x, n04 25x; n05 (10x) and n06 (25x) have another
         nuclear genome whose eight copies are all 99.9% identical, with substitutions only (recorded in truth).
         Organelle SNPs: n01 none; n02 = n05: 3; n03 = n02 + 2; n04 5; n06 = n04 + 2.

Reads mimic R10.4.1 SUP reads: ~0.4% substitutions, ~0.2% random indels and homopolymer length errors that grow
with the run length; read lengths are log-normal (N50 ~8 kb); 15% of the reads come from an unrelated sequence
(baiting must remove them). Deterministic: the same seed gives the same files. Standard library only. The
read names of heteroplasmy and numt say where each read comes from (a last token after "_": major, minor,
offtarget; organelle, nupt<insert id>, nuclear), for score.py to count them.

Writes, for each scenario: reference.fasta, reads/<sample>.fastq.gz, truth/<sample>.fasta (true genomes) and
truth.json (SNP events, expected distances, circularity, repeats); heteroplasmy adds truth/minor/<sample>.fasta
(the minor haplotype) and the minor sites, numt truth/nuclear/<A|B>.fasta and the inserts.

Usage: python simulate.py OUTPUT_FOLDER [--scenario plastid|linear|rdna|heteroplasmy|numt|all|new] [--depth 40]
       [--seed N]
'all' is plastid, linear and rdna (the scenarios of every validation record), 'new' heteroplasmy and numt. Each
scenario has its own random stream (seeded from the seed and its name), so adding scenarios changes none of the
others. --depth applies to plastid, linear and rdna; heteroplasmy and numt set their own depths.
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
    extra: dict = field(default_factory=dict)  # Scenario-specific truth (heteroplasmy, numt)


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


def _plastome(rng: random.Random, lsc: int, ir: int, ssc: int) -> tuple[str, list[dict]]:
    """LSC, IRa, SSC, IRb (= reverse complement of IRa), as in `plastid`, and the two repeats."""
    lsc_seq, ir_seq, ssc_seq = random_sequence(rng, lsc), random_sequence(rng, ir), random_sequence(rng, ssc)
    irb_start = lsc + ir + ssc
    repeats = [{"name": "IRa", "start": lsc, "end": lsc + ir}, {"name": "IRb", "start": irb_start,
                                                                "end": irb_start + ir}]
    return lsc_seq + ir_seq + ssc_seq + revcomp(ir_seq), repeats


# Minor-haplotype fraction of each sample; h020s also has SNPs of its own (in both haplotypes).
HETEROPLASMY = {"h000": 0.0, "h002": 0.02, "h005": 0.05, "h010": 0.10, "h020": 0.20, "h030": 0.30, "h040": 0.40,
                "h050": 0.50, "h020s": 0.20}


def heteroplasmy(rng: random.Random) -> Scenario:
    lsc, ir, ssc = 40_000, 10_000, 10_000
    genome, repeats = _plastome(rng, lsc, ir, ssc)
    irb_start = lsc + ir + ssc
    sc = Scenario("heteroplasmy", [Molecule("plastid", genome, True)], repeats=repeats)
    lsc_range, ssc_range = (300, lsc - 300), (lsc + ir + 300, irb_start - 300)
    used: set[int] = set()

    def snps(prefix: str, n: int, allowed: list[tuple[int, int]]) -> list[str]:
        ids = []
        for i, p in enumerate(_choose_positions(rng, n, allowed, used)):
            vid = f"{prefix}{i + 1}"
            sc.variants[vid] = Variant(vid, "snp", 0, [p], _alt(rng, genome[p]))
            ids.append(vid)
        return ids

    # The minor haplotype: 15 SNPs in the LSC, 5 in the SSC, 5 in IRa (mirrored in IRb) and three indels.
    minor = snps("mL", 15, [lsc_range]) + snps("mS", 5, [ssc_range])
    for i, off in enumerate(_choose_positions(rng, 5, [(lsc + 500, lsc + ir - 500)], used)):
        vid = f"mIR{i + 1}"
        mirror = irb_start + (ir - 1 - (off - lsc))
        sc.variants[vid] = Variant(vid, "snp", 0, [off, mirror], _alt(rng, genome[off]), in_repeat=True)
        minor.append(vid)
    for vid, kind, size in (("mdel1", "del", 1), ("mins1", "ins", 3), ("mdel2", "del", 5)):
        p = _choose_positions(rng, 1, [lsc_range, ssc_range], used)[0]
        if kind == "ins":
            sc.variants[vid] = Variant(vid, kind, 0, [p], random_sequence(rng, size))
        else:
            sc.variants[vid] = Variant(vid, kind, 0, [p], length=size)
        minor.append(vid)
    own = snps("o", 12, [lsc_range, ssc_range])
    sc.samples = {s: own if s.endswith("s") else [] for s in HETEROPLASMY}
    sc.depth = {s: 80.0 for s in HETEROPLASMY}
    sc.extra = {"minor_fraction": dict(HETEROPLASMY), "minor_variants": minor}
    return sc


# NUPTs: (length, identity to the plastome) of each insert, numbered from 1 in this order.
NUPT_INSERTS = ((500, 0.999), (500, 0.90), (2_000, 0.99), (2_000, 0.97), (5_000, 0.97), (5_000, 0.999),
                (15_000, 0.90), (15_000, 0.99))
NUPT_SPACER = 80_000  # Two of the longest reads (40 kb): no chain of baited reads joins two inserts
# Sample: (nuclear depth, nuclear genome). B: the same inserts, all 99.9% identical, substitutions only.
NUMT_SAMPLES = {"n01": (0.5, "A"), "n02": (2.0, "A"), "n03": (10.0, "A"), "n04": (25.0, "A"), "n05": (10.0, "B"),
                "n06": (25.0, "B")}
ORGANELLE_DEPTH = 50.0


def diverge(rng: random.Random, seq: str, identity: float, indels: bool) -> tuple[str, list[tuple[int, str, str]],
                                                                                    int]:
    """A copy of `seq` with round(length x (1 - identity)) changes (at least one): substitutions, and, if
    `indels`, 15% of 1-3 bp insertions or deletions. Returns the copy, its substitutions (offset, original base,
    new base) and the number of indels."""
    n = max(1, round(len(seq) * (1 - identity)))
    offsets = sorted(rng.sample(range(10, len(seq) - 10), n))
    out = list(seq)
    subs, n_indels = [], 0
    for i in range(n - 1, -1, -1):  # Right to left: offsets still point into the original sequence
        o = offsets[i]
        room = (offsets[i + 1] if i + 1 < n else len(seq)) - o - 1  # A deletion must not reach the next change
        if indels and rng.random() < 0.15 and room >= 1:
            size = rng.randint(1, 3)
            if rng.random() < 0.5:
                out[o + 1: o + 1] = list(random_sequence(rng, size, gc=0.40))
            else:
                del out[o: o + min(size, room)]
            n_indels += 1
        else:
            out[o] = _alt(rng, seq[o])
            subs.append((o, seq[o], out[o]))
    return "".join(out), sorted(subs), n_indels


def numt(rng: random.Random) -> Scenario:
    lsc, ir, ssc = 80_000, 12_000, 16_000
    genome, repeats = _plastome(rng, lsc, ir, ssc)
    irb_start = lsc + ir + ssc
    sc = Scenario("numt", [Molecule("plastid", genome, True)], repeats=repeats)
    single_copy = [(0, lsc), (lsc + ir, irb_start)]
    # Segments copied into the nucleus: in the single-copy regions, 1 kb apart; the longest placed first.
    sources: dict[int, tuple[int, int]] = {}
    for k in sorted(range(len(NUPT_INSERTS)), key=lambda i: -NUPT_INSERTS[i][0]):
        size = NUPT_INSERTS[k][0]
        regions = [(a, b) for a, b in single_copy if b - a - 1_000 > size]
        while True:
            a, b = rng.choice(regions)
            start = rng.randrange(a + 500, b - 500 - size)
            if all(start + size + 1_000 <= s or start >= e + 1_000 for s, e in sources.values()):
                sources[k] = (start, start + size)
                break
    # Organelle SNPs of the samples, anywhere in the single-copy regions (inside or outside the copied segments).
    used: set[int] = set()
    allowed = [(300, lsc - 300), (lsc + ir + 300, irb_start - 300)]

    def snps(prefix: str, n: int) -> list[str]:
        ids = []
        for i, p in enumerate(_choose_positions(rng, n, allowed, used)):
            vid = f"{prefix}{i + 1}"
            sc.variants[vid] = Variant(vid, "snp", 0, [p], _alt(rng, genome[p]))
            ids.append(vid)
        return ids

    a, b, c, d = snps("a", 3), snps("b", 2), snps("c", 5), snps("d", 2)
    sc.samples = {"n01": [], "n02": a, "n03": a + b, "n04": c, "n05": a, "n06": c + d}
    sc.depth = {s: ORGANELLE_DEPTH for s in sc.samples}
    # Nuclear genomes: random spacers and the inserts, on either strand; A and B share spacers, strands and
    # copied segments.
    spacers = [random_sequence(rng, NUPT_SPACER, gc=0.40) for _ in range(len(NUPT_INSERTS) + 1)]
    strands = ["-" if rng.random() < 0.5 else "+" for _ in NUPT_INSERTS]
    nuclear, inserts = {}, {}
    for name, fixed in (("A", None), ("B", 0.999)):
        parts, inserts[name] = [spacers[0]], []
        pos = NUPT_SPACER
        for k, (size, identity) in enumerate(NUPT_INSERTS):
            s, e = sources[k]
            copy, subs, n_indels = diverge(rng, genome[s:e], fixed or identity, indels=fixed is None)
            if strands[k] == "-":
                copy = revcomp(copy)
            inserts[name].append({"id": k + 1, "length": size, "identity": fixed or identity,
                                  "source_start": s, "source_end": e, "strand": strands[k],
                                  "nuclear_start": pos, "nuclear_end": pos + len(copy), "indels": n_indels,
                                  # Plastome position, plastome base, base of the copy (plastome strand)
                                  "substitutions": [[s + o, x, y] for o, x, y in subs]})
            parts += [copy, spacers[k + 1]]
            pos += len(copy) + NUPT_SPACER
        nuclear[name] = Molecule(f"nuclear_{name}", "".join(parts), False)
    sc.extra = {"nuclear": nuclear, "inserts": inserts, "copied_regions": sorted(sources.values()),
                "nuclear_depth": {s: dn for s, (dn, _) in NUMT_SAMPLES.items()},
                "nuclear_genome": {s: g for s, (_, g) in NUMT_SAMPLES.items()}}
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


def _fragment(mol: Molecule, length: int, rng: random.Random) -> tuple[int, str]:
    """A fragment of `mol` and its start (as in simulate_reads: circular molecules wrap, linear ones stop)."""
    if mol.circular:
        start = rng.randrange(len(mol.seq))
        return start, (mol.seq + mol.seq)[start: start + min(length, len(mol.seq))]
    start = rng.randrange(-length + 200, len(mol.seq) - 200)
    return max(0, start), mol.seq[max(0, start): start + length]


def _fastq(fragment: str, rng: random.Random, name: str) -> str:
    if rng.random() < 0.5:
        fragment = revcomp(fragment)
    seq = add_errors(fragment, rng)
    return f"@{name} simulated\n{seq}\n+\n{quality(rng, len(seq))}\n"


def mixture_reads(major: Molecule, minor: Molecule, fraction: float, off_target: str, depth: float,
                  rng: random.Random, name: str, off_target_fraction: float = 0.15) -> list[str]:
    """Reads of a mixture of two haplotypes, the minor one making `fraction` of the organelle reads, plus
    unrelated reads. Each read name ends in its origin: major, minor or offtarget."""
    records, bases, n = [], 0, 0
    while bases < depth * len(major.seq):
        n += 1
        length = read_length(rng)
        if rng.random() < off_target_fraction:
            start = rng.randrange(len(off_target) - length)
            fragment, origin = off_target[start: start + length], "offtarget"
        else:
            origin = "minor" if rng.random() < fraction else "major"
            fragment = _fragment(minor if origin == "minor" else major, length, rng)[1]
            if len(fragment) < 200:
                continue
            bases += len(fragment)
        records.append(_fastq(fragment, rng, f"{name}_{n}_{origin}"))
    return records


def numt_reads(organelle: Molecule, nuclear: Molecule, inserts: list[dict], organelle_depth: float,
               nuclear_depth: float, rng: random.Random, name: str) -> list[str]:
    """Organelle reads and nuclear reads (no other off-target reads), shuffled. Each read name ends in its
    origin: organelle, nupt<id> (a nuclear read covering at least 100 bp of insert <id>; the insert it covers
    most) or nuclear."""
    records, n = [], 0
    for mol, depth in ((organelle, organelle_depth), (nuclear, nuclear_depth)):
        bases = 0
        while bases < depth * len(mol.seq):
            start, fragment = _fragment(mol, read_length(rng), rng)
            if len(fragment) < 200:
                continue
            bases += len(fragment)
            n += 1
            origin = "organelle" if mol is organelle else "nuclear"
            if mol is nuclear:
                end = start + len(fragment)
                overlap, best = max((min(end, i["nuclear_end"]) - max(start, i["nuclear_start"]), i["id"])
                                    for i in inserts)
                if overlap >= 100:
                    origin = f"nupt{best}"
            records.append(_fastq(fragment, rng, f"{name}_{n}_{origin}"))
    rng.shuffle(records)
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


def write_reads(path: Path, records: list[str]) -> None:
    # mtime=0: the same seed gives byte-identical files.
    with gzip.GzipFile(path, "wb", compresslevel=6, mtime=0) as fh:
        fh.write("".join(records).encode())


def truth_record(sc: Scenario, depth: float) -> dict:
    return {
        "scenario": sc.name,
        "molecules": [{"name": m.name, "length": len(m.seq), "circular": m.circular} for m in sc.reference],
        "repeats": sc.repeats,
        "samples": {s: {"variants": v, "depth": sc.depth.get(s, depth)} for s, v in sc.samples.items()},
        "variants": {k: v.__dict__ for k, v in sc.variants.items()},
        "distances": expected_distances(sc),
    }


def build(sc: Scenario, out: Path, depth: float, rng: random.Random) -> None:
    (out / "reads").mkdir(parents=True, exist_ok=True)
    (out / "truth").mkdir(exist_ok=True)
    write_fasta(out / "reference.fasta", sc.reference)
    off_target = random_sequence(rng, 300_000, gc=0.45)
    for sample, variant_ids in sc.samples.items():
        genome = apply_variants(sc, variant_ids)
        write_fasta(out / "truth" / f"{sample}.fasta", genome)
        reads = simulate_reads(genome, off_target, sc.depth.get(sample, depth), rng, sample)
        write_reads(out / "reads" / f"{sample}.fastq.gz", reads)
    (out / "truth.json").write_text(json.dumps(truth_record(sc, depth), indent=1) + "\n")


def minor_sites(sc: Scenario, major: str) -> list[dict]:
    """Every position where the minor haplotype differs from the major one, in major (= reference) coordinates:
    a SNP in the inverted repeat gives two sites. ins: `minor` is inserted after `position`; del: `major` (the
    bases from `position`) is deleted."""
    sites = []
    for vid in sc.extra["minor_variants"]:
        v = sc.variants[vid]
        if v.kind == "snp":
            for k, p in enumerate(v.positions):
                base = v.alt if k == 0 else v.alt.translate(COMPLEMENT)
                sites.append({"id": vid if k == 0 else f"{vid}b", "kind": "snp", "position": p,
                              "major": major[p], "minor": base, "in_repeat": v.in_repeat})
        elif v.kind == "ins":
            sites.append({"id": vid, "kind": "ins", "position": v.positions[0], "major": "", "minor": v.alt,
                          "in_repeat": False})
        else:
            p = v.positions[0]
            sites.append({"id": vid, "kind": "del", "position": p, "major": major[p: p + v.length], "minor": "",
                          "in_repeat": False})
    return sites


def build_heteroplasmy(sc: Scenario, out: Path, rng: random.Random) -> None:
    (out / "reads").mkdir(parents=True, exist_ok=True)
    (out / "truth" / "minor").mkdir(parents=True, exist_ok=True)
    write_fasta(out / "reference.fasta", sc.reference)
    off_target = random_sequence(rng, 300_000, gc=0.45)
    truth = truth_record(sc, 0.0)
    for sample, variant_ids in sc.samples.items():
        major = apply_variants(sc, variant_ids)
        minor = apply_variants(sc, variant_ids + sc.extra["minor_variants"])
        write_fasta(out / "truth" / f"{sample}.fasta", major)
        write_fasta(out / "truth" / "minor" / f"{sample}.fasta", minor)
        fraction = sc.extra["minor_fraction"][sample]
        reads = mixture_reads(major[0], minor[0], fraction, off_target, sc.depth[sample], rng, sample)
        write_reads(out / "reads" / f"{sample}.fastq.gz", reads)
        truth["samples"][sample].update(minor_fraction=fraction, minor_sites=minor_sites(sc, major[0].seq))
    truth["minor_variants"] = sc.extra["minor_variants"]
    (out / "truth.json").write_text(json.dumps(truth, indent=1) + "\n")


def build_numt(sc: Scenario, out: Path, rng: random.Random) -> None:
    (out / "reads").mkdir(parents=True, exist_ok=True)
    (out / "truth" / "nuclear").mkdir(parents=True, exist_ok=True)
    write_fasta(out / "reference.fasta", sc.reference)
    for name, mol in sc.extra["nuclear"].items():
        write_fasta(out / "truth" / "nuclear" / f"{name}.fasta", [mol])
    truth = truth_record(sc, 0.0)
    for sample, variant_ids in sc.samples.items():
        organelle = apply_variants(sc, variant_ids)
        write_fasta(out / "truth" / f"{sample}.fasta", organelle)
        genome = sc.extra["nuclear_genome"][sample]
        depth = sc.extra["nuclear_depth"][sample]
        reads = numt_reads(organelle[0], sc.extra["nuclear"][genome], sc.extra["inserts"][genome],
                           sc.depth[sample], depth, rng, sample)
        write_reads(out / "reads" / f"{sample}.fastq.gz", reads)
        truth["samples"][sample].update(nuclear_depth=depth, nuclear_genome=genome)
    truth["copied_regions"] = sc.extra["copied_regions"]
    truth["nuclear"] = {name: {"length": len(mol.seq), "inserts": sc.extra["inserts"][name]}
                        for name, mol in sc.extra["nuclear"].items()}
    (out / "truth.json").write_text(json.dumps(truth, indent=1) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("output", type=Path)
    parser.add_argument("--scenario", choices=["plastid", "linear", "rdna", "heteroplasmy", "numt", "all", "new"],
                        default="all", help="all: plastid, linear and rdna; new: heteroplasmy and numt")
    parser.add_argument("--depth", type=float, default=40.0, help="plastid, linear and rdna only")
    parser.add_argument("--seed", type=int, default=20260924)
    args = parser.parse_args()
    for name, make in (("plastid", plastid), ("linear", linear), ("rdna", rdna)):
        if args.scenario in (name, "all"):
            rng = random.Random(f"{args.seed}-{name}")
            build(make(rng), args.output / name, args.depth, rng)
            print(f"{name}: written to {args.output / name}")
    for name, make, build_new in (("heteroplasmy", heteroplasmy, build_heteroplasmy), ("numt", numt, build_numt)):
        if args.scenario in (name, "new"):
            rng = random.Random(f"{args.seed}-{name}")
            build_new(make(rng), args.output / name, rng)
            print(f"{name}: written to {args.output / name}")


if __name__ == "__main__":
    main()
