#!/usr/bin/env python3
"""Generate the BACoN example: an annotated 30 kb circular "organelle" reference, four simulated nanopore samples
with known SNPs, and a sample metadata file, so that every part of the report is exercised.

The reference is laid out like a plastome: a large single-copy region (LSC, 16,000 bp), an inverted repeat
(IRb, 3,000 bp), a small single-copy region (SSC, 8,000 bp) and the second copy of the repeat (IRa, the reverse
complement of IRb). It carries 15 synthetic protein-coding genes (orf01..orf15: open reading frames of random
sense codons on both strands, orf04 split by an intron), tRNA genes, an rRNA gene and a tRNA gene in each copy of
the inverted repeat, and a pseudogene (a truncated, diverged copy of orf03). Every name is made up: nothing in it
comes from a real organism. reference.gb is its GenBank record and reference.fasta the same sequence as fasta,
byte for byte what BACoN writes from the GenBank file.

Each sample carries known SNPs relative to the reference, so the expected SNP distances are known:

    sample   SNPs
    alpha    none (identical to the reference)
    beta     6 SNPs
    gamma    beta's 6 SNPs + 4 of its own
    delta    10 SNPs of its own

The SNPs are planted outside the inverted repeats (both copies stay identical in every genome), at least 40 bp
apart, in chosen contexts: coding sequences (synonymous, missense and nonsense changes, on both strands, one in
the second exon of the split gene), an intron, tRNAs, the pseudogene and intergenic spacers. planted_snps.tsv
lists them per sample; planted_effects.tsv gives the expected region, gene, context and effect of each, computed
here without BACoN. metadata.tsv is a sample metadata file with simulated values.

Reads are 1-20 kb long with about 1.5% errors (substitutions, insertions, deletions), at ~50x depth, mixed with
20% off-target reads from an unrelated sequence, which baiting must remove. Deterministic (fixed seed).
Standard library only.

Usage: python make_example.py OUTPUT_FOLDER
"""

from __future__ import annotations

import gzip
import random
import sys
from dataclasses import dataclass, field
from pathlib import Path

GENOME_LENGTH = 30_000
DEPTH = 50
ERROR_RATE = 0.015
OFF_TARGET_FRACTION = 0.2
SEED = 20260924
MIN_SNP_SPACING = 40  # SKA's split k-mers need clean flanks
MIN_END_DISTANCE = 500

# The regions (1-based, inclusive): LSC, IRb, SSC, IRa; IRa is the reverse complement of IRb.
LSC_LENGTH, IR_LENGTH, SSC_LENGTH = 16_000, 3_000, 8_000
REGIONS = {"LSC": (1, LSC_LENGTH), "IRb": (LSC_LENGTH + 1, LSC_LENGTH + IR_LENGTH),
           "SSC": (LSC_LENGTH + IR_LENGTH + 1, LSC_LENGTH + IR_LENGTH + SSC_LENGTH),
           "IRa": (LSC_LENGTH + IR_LENGTH + SSC_LENGTH + 1, GENOME_LENGTH)}
assert LSC_LENGTH + 2 * IR_LENGTH + SSC_LENGTH == GENOME_LENGTH

SAMPLES = {"alpha": [], "beta": ["b"], "gamma": ["b", "g"], "delta": ["d"]}  # The SNP sets of each sample

EXPECTED_DISTANCES = {
    ("alpha", "beta"): 6, ("alpha", "gamma"): 10, ("alpha", "delta"): 10,
    ("beta", "gamma"): 4, ("beta", "delta"): 16, ("gamma", "delta"): 20,
}

# The genetic code (table 11: the standard code with the bacterial and plastid start codons; only ATG is used).
_CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = {a + b + c: _CODE[16 * i + 4 * j + k]
          for i, a in enumerate("TCAG") for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")}
SENSE_CODONS = sorted(c for c, aa in CODONS.items() if aa != "*")
STOP_CODONS = sorted(c for c, aa in CODONS.items() if aa == "*")


@dataclass
class Feature:
    """An annotated gene (CDS, tRNA, rRNA or pseudogene): its parts in coordinate order, 1-based, inclusive."""

    name: str
    kind: str
    strand: int
    parts: list[tuple[int, int]]
    product: str
    note: str = ""
    translation: str = ""

    @property
    def start(self) -> int:
        return self.parts[0][0]

    @property
    def end(self) -> int:
        return self.parts[-1][1]

    def contains(self, pos: int) -> bool:
        return self.start <= pos <= self.end

    def in_exon(self, pos: int) -> bool:
        return any(s <= pos <= e for s, e in self.parts)

    def coding_index(self, pos: int) -> int:
        """0-based index of a reference position in the coding sequence (5' to 3' on the coding strand)."""
        ordered = self.parts[::-1] if self.strand < 0 else self.parts
        offset = 0
        for s, e in ordered:
            if s <= pos <= e:
                return offset + (e - pos if self.strand < 0 else pos - s)
            offset += e - s + 1
        raise ValueError(pos)

    def position_of(self, index: int) -> int:
        """The reference position of a coding-sequence index (the inverse of coding_index)."""
        ordered = self.parts[::-1] if self.strand < 0 else self.parts
        for s, e in ordered:
            if index < e - s + 1:
                return e - index if self.strand < 0 else s + index
            index -= e - s + 1
        raise ValueError(index)

    def location(self, hull: bool = False) -> str:
        """The GenBank location of the parts, or of their hull (the gene feature of a split gene)."""
        parts = [(self.start, self.end)] if hull else self.parts
        text = ",".join(f"{s}..{e}" for s, e in parts)
        if len(parts) > 1:
            text = f"join({text})"
        return f"complement({text})" if self.strand < 0 else text


# The genes of the single-copy regions: name, kind, strand, parts. Lengths of the coding sequences are multiples
# of 3 (orf04: 400 + 599 bp). The gaps between them are intergenic spacers.
GENE_LAYOUT = [
    ("orf01", "CDS", 1, [(301, 1200)]),
    ("trnA-sim", "tRNA", 1, [(1401, 1475)]),
    ("orf02", "CDS", -1, [(1601, 2800)]),
    ("orf03", "CDS", 1, [(3001, 3900)]),
    ("orf04", "CDS", 1, [(4201, 4600), (5101, 5699)]),  # Split by an intron (4601..5100)
    ("orf05", "CDS", -1, [(5901, 6500)]),
    ("orf06", "CDS", 1, [(6701, 8200)]),
    ("orf03-ps", "pseudogene", 1, [(8401, 8700)]),  # Truncated, diverged copy of the start of orf03
    ("orf07", "CDS", -1, [(8901, 9800)]),
    ("trnB-sim", "tRNA", -1, [(10001, 10075)]),
    ("orf08", "CDS", 1, [(10201, 11400)]),
    ("orf09", "CDS", 1, [(11601, 11900)]),
    ("orf10", "CDS", -1, [(12101, 13300)]),
    ("orf11", "CDS", 1, [(13501, 14400)]),
    ("orf12", "CDS", 1, [(19301, 20200)]),
    ("orf13", "CDS", -1, [(20401, 21900)]),
    ("trnC-sim", "tRNA", 1, [(22101, 22175)]),
    ("orf14", "CDS", 1, [(22401, 23000)]),
    ("orf15", "CDS", -1, [(24001, 24900)]),
]
# The genes of the inverted repeat, in IRb; IRa carries their mirror images on the other strand.
IR_GENE_LAYOUT = [
    ("rrn16-sim", "rRNA", 1, [(16301, 17800)]),
    ("trnD-sim", "tRNA", 1, [(18001, 18075)]),
]
PSEUDOGENE_SOURCE = "orf03"

# The planted SNPs: SNP set, where, and what. In a CDS, the first position from the given fraction of the coding
# sequence whose change gives the intended effect is used (missense and nonsense changes keep the start codon);
# an intron, tRNA or pseudogene SNP is at the given fraction of the feature; an intergenic SNP at the position.
SNP_PLAN = [
    ("b", "orf01", "synonymous", 0.3),
    ("b", "orf02", "missense", 0.3),  # Reverse strand
    ("b", "orf04", "missense", 0.7),  # Second exon of the split gene
    ("b", "orf04", "intron", 0.5),
    ("b", "intergenic", "LSC", 15_000),
    ("b", "orf13", "synonymous", 0.5),  # Reverse strand, SSC
    ("g", "orf06", "missense", 0.4),
    ("g", "trnB-sim", "tRNA", 0.5),
    ("g", "intergenic", "SSC", 23_500),
    ("g", "orf12", "synonymous", 0.6),
    ("d", "orf01", "missense", 0.7),
    ("d", "orf03", "nonsense", 0.5),
    ("d", "orf05", "synonymous", 0.4),  # Reverse strand
    ("d", "orf07", "missense", 0.6),  # Reverse strand
    ("d", "orf03-ps", "pseudogene", 0.5),
    ("d", "intergenic", "LSC", 12_000),
    ("d", "trnC-sim", "tRNA", 0.5),
    ("d", "orf14", "missense", 0.3),
    ("d", "orf15", "synonymous", 0.5),  # Reverse strand, SSC
    ("d", "intergenic", "SSC", 26_000),
]


@dataclass
class PlantedSnp:
    position: int  # 1-based
    ref: str
    alt: str
    sets: list[str]
    region: str
    gene: str
    context: str
    codon: str = ""  # REF>ALT codon on the coding strand
    aa: str = ""  # Amino-acid change (M1T)
    effect: str = ""

    @property
    def samples(self) -> list[str]:
        return [s for s, sets in SAMPLES.items() if any(k in sets for k in self.sets)]


@dataclass
class Reference:
    seq: str
    features: list[Feature]
    snps: list[PlantedSnp] = field(default_factory=list)

    def snp_sets(self) -> dict[str, list[tuple[int, str]]]:
        """The (position, alt) of each SNP set."""
        sets: dict[str, list[tuple[int, str]]] = {}
        for snp in self.snps:
            for k in snp.sets:
                sets.setdefault(k, []).append((snp.position, snp.alt))
        return {k: sorted(v) for k, v in sets.items()}


# ---------------------------------------------------------------------------------------------------------------
# Sequences
# ---------------------------------------------------------------------------------------------------------------

def random_sequence(rng: random.Random, length: int, gc: float = 0.37) -> str:
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights=weights, k=length))


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def translate(coding: str) -> str:
    return "".join(CODONS[coding[i:i + 3]] for i in range(0, len(coding) - 2, 3))


def make_orf(rng: random.Random, length: int) -> str:
    """An open reading frame of `length` bp: ATG, random sense codons, a stop codon."""
    assert length % 3 == 0 and length >= 9
    return "ATG" + "".join(rng.choices(SENSE_CODONS, k=length // 3 - 2)) + rng.choice(STOP_CODONS)


def coding_sequence(seq: str, feature: Feature) -> str:
    text = "".join(seq[s - 1:e] for s, e in feature.parts)
    return revcomp(text) if feature.strand < 0 else text


def _place(genome: list[str], feature: Feature, coding: str) -> None:
    """Write a coding sequence into the genome over the parts of its feature."""
    text = revcomp(coding) if feature.strand < 0 else coding
    assert len(text) == sum(e - s + 1 for s, e in feature.parts)
    i = 0
    for s, e in feature.parts:
        genome[s - 1:e] = text[i:i + e - s + 1]
        i += e - s + 1


def region_of(pos: int) -> str:
    return next(name for name, (s, e) in REGIONS.items() if s <= pos <= e)


def make_reference(rng: random.Random) -> Reference:
    """The reference: random background, the genes written over it, the inverted repeat mirrored."""
    genome = list(random_sequence(rng, GENOME_LENGTH))
    features: list[Feature] = []
    orfs: dict[str, str] = {}
    number = 0
    for name, kind, strand, parts in GENE_LAYOUT:
        if kind == "CDS":
            number += 1
            orf = make_orf(rng, sum(e - s + 1 for s, e in parts))
            feature = Feature(name, kind, strand, parts, f"simulated protein {number}",
                              translation=translate(orf)[:-1])
            _place(genome, feature, orf)
            orfs[name] = orf
        elif kind == "pseudogene":
            # The first bases of the source gene with every sixth base changed (so that no 31-mer is shared with
            # the source) and an early stop codon
            length = sum(e - s + 1 for s, e in parts)
            copy = list(orfs[PSEUDOGENE_SOURCE][:length])
            for i in range(5, length, 6):
                copy[i] = "ACGT"[("ACGT".index(copy[i]) + 1) % 4]
            copy[60:63] = "TAA"
            feature = Feature(name, kind, strand, parts, f"pseudogene, truncated copy of {PSEUDOGENE_SOURCE}",
                              note=f"truncated and diverged copy of the 5' end of {PSEUDOGENE_SOURCE}")
            _place(genome, feature, "".join(copy))
        else:
            feature = Feature(name, kind, strand, parts, f"simulated {kind} {name.split('-')[0][3:]}")
        features.append(feature)
    irb_start, irb_end = REGIONS["IRb"]
    ira_start, ira_end = REGIONS["IRa"]
    for name, kind, strand, parts in IR_GENE_LAYOUT:
        product = f"simulated {kind} {name.split('-')[0][3:]}"
        features.append(Feature(name, kind, strand, parts, product, note="in IRb"))
        mirrored = sorted((ira_start + irb_end - e, ira_start + irb_end - s) for s, e in parts)
        features.append(Feature(name, kind, -strand, mirrored, product, note="in IRa"))
    genome[ira_start - 1:ira_end] = revcomp("".join(genome[irb_start - 1:irb_end]))
    features.sort(key=lambda f: (f.start, f.end))
    reference = Reference("".join(genome), features)
    reference.snps = plant_snps(reference, rng)
    return reference


# ---------------------------------------------------------------------------------------------------------------
# Planted SNPs and their expected effects
# ---------------------------------------------------------------------------------------------------------------

def _cds_change(coding: str, index: int, alt: str) -> tuple[str, str, str, str]:
    """The codon change, amino-acid change and effect of a base change in a coding sequence."""
    codon_index, within = divmod(index, 3)
    codon = coding[3 * codon_index:3 * codon_index + 3]
    codon_alt = codon[:within] + alt + codon[within + 1:]
    aa, aa_alt = CODONS[codon], CODONS[codon_alt]
    if codon_index == 0:
        effect = "start retained" if codon_alt == "ATG" else "start lost"
    elif aa == aa_alt:
        effect = "stop retained" if aa == "*" else "synonymous"
    elif aa == "*":
        effect = "stop lost"
    elif aa_alt == "*":
        effect = "nonsense"
    else:
        effect = "missense"
    return f"{codon}>{codon_alt}", f"{aa}{codon_index + 1}{aa_alt}", effect, codon_alt


def _cds_snp(seq: str, feature: Feature, effect: str, fraction: float) -> tuple[int, str, str, str]:
    """The first change from `fraction` of the coding sequence with the intended effect: position, alternate base
    in reference coordinates, codon change, amino-acid change."""
    coding = coding_sequence(seq, feature)
    for index in range(max(3, int(len(coding) * fraction)), len(coding) - 3):
        for alt in "ACGT":
            if alt == coding[index]:
                continue
            codons, aa, kind, _ = _cds_change(coding, index, alt)
            if kind == effect:
                base = revcomp(alt) if feature.strand < 0 else alt
                return feature.position_of(index), base, codons, aa
    raise ValueError(f"no {effect} change in {feature.name} from {fraction}")


def _neighbours(features: list[Feature], pos: int) -> tuple[Feature, Feature]:
    """The nearest gene ending before and the nearest starting after an intergenic position (circular)."""
    before = max((f for f in features if f.end < pos), key=lambda f: f.end, default=None) or \
        max(features, key=lambda f: f.end)
    after = min((f for f in features if f.start > pos), key=lambda f: f.start, default=None) or \
        min(features, key=lambda f: f.start)
    return before, after


def plant_snps(reference: Reference, rng: random.Random) -> list[PlantedSnp]:
    seq, features = reference.seq, reference.features
    by_name = {f.name: f for f in features if f.kind != "rRNA"}
    planted: dict[int, PlantedSnp] = {}
    for key, where, what, value in SNP_PLAN:
        codon = aa = effect = ""
        if where == "intergenic":
            pos = int(value)
            before, after = _neighbours(features, pos)
            gene, context = "", f"intergenic between {before.name} and {after.name}"
            assert not any(f.contains(pos) for f in features), pos
        else:
            feature = by_name[where]
            gene = feature.name
            if what in ("synonymous", "missense", "nonsense"):
                pos, alt, codon, aa = _cds_snp(seq, feature, what, value)
                context, effect = "CDS", what
            elif what == "intron":
                pos = feature.parts[0][1] + 1 + int((feature.parts[1][0] - feature.parts[0][1] - 1) * value)
                context = "intron"
            else:
                pos = feature.start + int((feature.end - feature.start) * value)
                context = what
        if not effect:
            alt = rng.choice([b for b in "ACGT" if b != seq[pos - 1]])
        if pos in planted:
            planted[pos].sets.append(key)
            continue
        planted[pos] = PlantedSnp(pos, seq[pos - 1], alt, [key], region_of(pos), gene, context, codon, aa, effect)
    snps = sorted(planted.values(), key=lambda s: s.position)
    positions = [s.position for s in snps]
    assert all(b - a >= MIN_SNP_SPACING for a, b in zip(positions, positions[1:])), positions
    assert positions[0] >= MIN_END_DISTANCE and positions[-1] <= GENOME_LENGTH - MIN_END_DISTANCE, positions
    assert all(s.region in ("LSC", "SSC") for s in snps), positions
    return snps


def mutate(seq: str, changes: list[tuple[int, str]]) -> str:
    s = list(seq)
    for pos, alt in changes:
        s[pos - 1] = alt
    return "".join(s)


# ---------------------------------------------------------------------------------------------------------------
# Output files
# ---------------------------------------------------------------------------------------------------------------

def fasta_text(reference: Reference) -> str:
    lines = [">organelle simulated circular reference"]
    lines += [reference.seq[j:j + 80] for j in range(0, len(reference.seq), 80)]
    return "\n".join(lines) + "\n"


def _qualifier(key: str, value: str | None = None, quoted: bool = True) -> list[str]:
    """A GenBank qualifier, wrapped at 79 columns (21 of indent): unquoted, or quoted and wrapped at spaces, or a
    translation wrapped anywhere."""
    indent = " " * 21
    if value is None:
        return [f"{indent}/{key}"]
    if not quoted:
        return [f"{indent}/{key}={value}"]
    width = 79 - len(indent)
    text = f'/{key}="{value}"'
    if key == "translation":
        return [indent + text[i:i + width] for i in range(0, len(text), width)]
    lines, current = [], ""
    for word in text.split(" "):
        if current and len(current) + 1 + len(word) > width:
            lines.append(current)
            current = word
        else:
            current = f"{current} {word}" if current else word
    lines.append(current)
    return [indent + line for line in lines]


def genbank_text(reference: Reference) -> str:
    """The GenBank record of the reference: LOCUS, DEFINITION, VERSION, the features and the sequence."""
    seq = reference.seq
    lines = [f"LOCUS       {'organelle':<16} {len(seq):>11} bp    DNA     circular SYN 01-JAN-2026",
             "DEFINITION  simulated circular reference.",
             "ACCESSION   organelle",
             "VERSION     organelle",
             "KEYWORDS    .",
             "SOURCE      synthetic construct",
             "  ORGANISM  synthetic construct",
             "            other sequences; artificial sequences.",
             "FEATURES             Location/Qualifiers",
             f"     source          1..{len(seq)}"]
    lines += _qualifier("organism", "synthetic construct")
    lines += _qualifier("mol_type", "genomic DNA")
    lines += _qualifier("note", "simulated circular reference of the BACoN example, generated by "
                                "example/make_example.py; every gene is synthetic")
    for name, (s, e) in REGIONS.items():
        if name.startswith("IR"):
            lines.append(f"     repeat_region   {s}..{e}")
            lines += _qualifier("rpt_type", "inverted", quoted=False)
            lines += _qualifier("standard_name", name)
            lines += _qualifier("note", f"inverted repeat {name[-1].upper()}")
    for f in reference.features:
        lines.append(f"     gene            {f.location(hull=True)}")
        lines += _qualifier("gene", f.name)
        if f.kind == "pseudogene":
            lines += _qualifier("pseudo")
            lines += _qualifier("note", f.note)
            continue
        lines.append(f"     {f.kind:<16}{f.location()}")
        lines += _qualifier("gene", f.name)
        if f.kind == "CDS":
            lines += _qualifier("codon_start", "1", quoted=False)
            lines += _qualifier("transl_table", "11", quoted=False)
        lines += _qualifier("product", f.product)
        if f.note:
            lines += _qualifier("note", f.note)
        if f.kind == "CDS":
            lines += _qualifier("translation", f.translation)
    lines.append("ORIGIN")
    for i in range(0, len(seq), 60):
        lines.append(f"{i + 1:>9} " + " ".join(seq[j:j + 10].lower() for j in range(i, min(i + 60, len(seq)), 10)))
    lines.append("//")
    return "\n".join(lines) + "\n"


def planted_snps_text(reference: Reference) -> str:
    sets = reference.snp_sets()
    lines = ["sample\tpositions (1-based)"]
    for sample, keys in SAMPLES.items():
        lines.append(f"{sample}\t{','.join(str(p) for k in keys for p, _ in sets[k])}")
    return "\n".join(lines) + "\n"


def planted_effects_text(reference: Reference) -> str:
    lines = ["position\tref\talt\tsamples\tregion\tgene\tcontext\tcodon\tamino_acid\teffect"]
    for s in reference.snps:
        lines.append("\t".join([str(s.position), s.ref, s.alt, ",".join(s.samples), s.region, s.gene, s.context,
                                s.codon, s.aa, s.effect]))
    return "\n".join(lines) + "\n"


METADATA = """\
# Simulated sample metadata for the BACoN example: every value is invented. 'group' colours the figures (the first
# column that can); delta has no group (NA); 'note' is free text.
sample\tgroup\tyear\torigin\tnote
alpha\tA\t2021\tsite 1\tIdentical to the reference by design (no planted SNP)
beta\tA\t2021\tsite 1\tCarries six planted SNPs, in coding sequences, an intron and a spacer
gamma\tB\t2022\tsite 2\tShares beta's six SNPs, plus four of its own (one in a tRNA)
delta\tNA\t2023\tsite 2\tTen SNPs of its own, including a nonsense change and one in the pseudogene
"""


def write_reference(out: Path, reference: Reference) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "reference.fasta").write_text(fasta_text(reference))
    (out / "reference.gb").write_text(genbank_text(reference))
    (out / "planted_snps.tsv").write_text(planted_snps_text(reference))
    (out / "planted_effects.tsv").write_text(planted_effects_text(reference))
    (out / "metadata.tsv").write_text(METADATA)


# ---------------------------------------------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------------------------------------------

def add_errors(seq: str, rng: random.Random) -> str:
    out = []
    for base in seq:
        r = rng.random()
        if r < ERROR_RATE * 0.5:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
        elif r < ERROR_RATE * 0.75:
            out.append(base + rng.choice("ACGT"))
        elif r < ERROR_RATE:
            continue
        else:
            out.append(base)
    return "".join(out)


def read_length(rng: random.Random) -> int:
    return max(1000, min(20_000, int(rng.lognormvariate(8.6, 0.5))))


def simulate(genome: str, off_target: str, rng: random.Random, name: str) -> list[str]:
    records, bases, n = [], 0, 0
    circular = genome + genome
    while bases < DEPTH * len(genome):
        n += 1
        length = read_length(rng)
        if rng.random() < OFF_TARGET_FRACTION:
            start = rng.randrange(len(off_target) - length)
            fragment = off_target[start: start + length]
        else:
            start = rng.randrange(len(genome))
            fragment = circular[start: start + min(length, len(genome))]
            bases += len(fragment)
        if rng.random() < 0.5:
            fragment = revcomp(fragment)
        seq = add_errors(fragment, rng)
        qual = "".join(rng.choice("0123456789:;<=>?@ABC") for _ in seq)
        records.append(f"@{name}_read{n} simulated\n{seq}\n+\n{qual}\n")
    return records


def write_reads(out: Path, reference: Reference, rng: random.Random) -> None:
    off_target = random_sequence(rng, 200_000, gc=0.5)
    sets = reference.snp_sets()
    (out / "reads").mkdir(parents=True, exist_ok=True)
    for sample, keys in SAMPLES.items():
        genome = mutate(reference.seq, [change for k in keys for change in sets[k]])
        with gzip.GzipFile(out / "reads" / f"{sample}.fastq.gz", "wb", compresslevel=6, mtime=0) as fh:
            fh.write("".join(simulate(genome, off_target, rng, sample)).encode())


def main(out: Path) -> None:
    rng = random.Random(SEED)
    reference = make_reference(rng)
    write_reference(out, reference)
    write_reads(out, reference, rng)
    print(f"Example written to {out}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(Path(sys.argv[1]))
