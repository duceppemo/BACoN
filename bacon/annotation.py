"""Reference annotations (GenBank or GFF3): the genes and regions of the reference, and the effect of the VCF's
SNPs on its coding sequences. Used by the report only (standard library only).

Sequence names: a GenBank record is named after its VERSION (accession.version, `NC_008096.2`), or its LOCUS name
without one, as NCBI's fasta of the same record; a GFF3 feature after its first column. Translation: each CDS's
`/transl_table` (GFF3: `transl_table=`), else table 11 (bacterial and plastid code); the genetic codes 1, 4 and 11
are known, other tables are translated with the standard code and table 11's start codons.
"""

from __future__ import annotations

import logging
import re
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote

from bacon import BaconError
from bacon.seqio import DECOMPRESSION_ERRORS, Record, open_text

log = logging.getLogger(__name__)

GENBANK_EXTENSIONS = (".gb", ".gbk", ".gbff", ".genbank")
GFF_EXTENSIONS = (".gff", ".gff3")
DEFAULT_TABLE = 11
MIN_INVERTED_REPEAT = 500  # Shorter inverted repeats (hairpins, transposon ends) do not make a region band
RNA_TYPES = {"tRNA": "tRNA", "rRNA": "rRNA", "ncRNA": "ncRNA", "misc_RNA": "ncRNA", "tmRNA": "ncRNA",
             "transcript": "ncRNA", "antisense_RNA": "ncRNA", "RNase_P_RNA": "ncRNA", "SRP_RNA": "ncRNA"}
REGION_TYPES = {"repeat_region", "misc_feature", "region", "inverted_repeat", "sequence_feature", "misc_structure",
                "biological_region", "repeat_unit"}
_BIN = 10_000  # Genes are indexed by bins of this many bases


def annotation_format(path: Path) -> str | None:
    """'genbank', 'gff3' or None, from the file name (.gb, .gbk, .gbff, .genbank, .gff, .gff3, with or without
    .gz), or else from its first line (LOCUS, ##gff-version)."""
    lower = path.name.lower()
    if lower.endswith(".gz"):
        lower = lower[:-3]
    if lower.endswith(GENBANK_EXTENSIONS):
        return "genbank"
    if lower.endswith(GFF_EXTENSIONS):
        return "gff3"
    try:
        with open_text(path) as fh:
            for line in fh:
                if line.strip():
                    if line.startswith("LOCUS "):
                        return "genbank"
                    if line.startswith("##gff-version"):
                        return "gff3"
                    return None
    except (OSError, *DECOMPRESSION_ERRORS):
        return None
    return None


# ---------------------------------------------------------------------------------------------------------------
# Locations
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class Location:
    strand: int  # 1 or -1
    parts: list[tuple[int, int]]  # 1-based, inclusive, in ascending order
    partial_low: bool = False  # '<' before the lowest coordinate
    partial_high: bool = False  # '>' after the highest coordinate
    listed: list[tuple[int, int, int]] | None = None  # GenBank: (start, end, strand) in the order written

    @property
    def ordered(self) -> list[tuple[int, int, int]]:
        """The parts, with their strand, in the order of translation: as the GenBank location lists them (a
        trans-spliced CDS has them out of coordinate order, or on both strands), else by coordinate, descending
        on the reverse strand."""
        if self.listed is not None:
            return self.listed
        return [(s, e, self.strand) for s, e in (self.parts[::-1] if self.strand < 0 else self.parts)]

    @property
    def start(self) -> int:
        return self.parts[0][0]

    @property
    def end(self) -> int:
        return self.parts[-1][1]


_RANGE = re.compile(r"^(<?)(\d+)(?:(?:\.\.|\^|\.)(>?)(\d+))?$")


def _split_top(text: str) -> list[str]:
    """Split on the commas outside parentheses."""
    items, depth, start = [], 0, 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "," and depth == 0:
            items.append(text[start:i])
            start = i + 1
    items.append(text[start:])
    return items


def _location_pieces(text: str, strand: int = 1) -> list[tuple[int, int, int, bool, bool]]:
    if text.startswith("complement(") and text.endswith(")"):
        return _location_pieces(text[11:-1], -strand)[::-1]  # Listed 5' to 3' of the complemented strand
    for prefix in ("join(", "order("):
        if text.startswith(prefix) and text.endswith(")"):
            pieces = []
            for item in _split_top(text[len(prefix):-1]):
                pieces += _location_pieces(item, strand)
            return pieces
    if ":" in text:  # A remote reference (ACCESSION:range): not on this sequence
        return []
    match = _RANGE.match(text)
    if not match:
        raise ValueError(text)
    start = int(match.group(2))
    end = int(match.group(4)) if match.group(4) else start
    if end < start or start < 1:
        raise ValueError(text)
    return [(start, end, strand, bool(match.group(1)), bool(match.group(3)))]


def parse_location(text: str) -> Location | None:
    """A GenBank location: `123..456`, `complement(...)`, `join(...)`, `order(...)`, nested, with `<`/`>` partial
    ends, single bases and `n^m` sites. None when every part is on another sequence (remote references, skipped).
    ValueError for anything else."""
    pieces = _location_pieces("".join(text.split()))
    if not pieces:
        return None
    listed = [(s, e, strand) for s, e, strand, *_ in pieces]
    pieces.sort()
    # The strand of the feature: that of most of its bases (a trans-spliced CDS can have parts on both strands)
    strand = 1 if sum((e - s + 1) * strand for s, e, strand in listed) >= 0 else -1
    return Location(strand, [(s, e) for s, e, *_ in pieces], pieces[0][3], pieces[-1][4], listed)


# ---------------------------------------------------------------------------------------------------------------
# GenBank and GFF3 reading
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class RawFeature:
    """A feature of either format, before features are grouped into genes."""

    type: str
    seq: str
    location: Location
    qualifiers: dict[str, str] = field(default_factory=dict)  # Lower-case keys; repeated values joined by '; '
    id: str = ""  # GFF3 ID
    parents: tuple[str, ...] = ()  # GFF3 Parent
    phases: dict[int, int] = field(default_factory=dict)  # GFF3: phase of the CDS part starting at this position


@dataclass
class GenBankRecord:
    name: str  # VERSION, or the LOCUS name
    definition: str
    length: int
    circular: bool
    features: list[RawFeature]
    seq: str
    skipped: int = 0  # Features whose location could not be parsed

    @property
    def header(self) -> str:
        """The fasta header NCBI gives the record: name and definition."""
        return f"{self.name} {self.definition}" if self.definition else self.name


def _add_qualifier(qualifiers: dict[str, str], key: str, value: str) -> None:
    key = key.lower()
    if key in qualifiers and value:
        qualifiers[key] = f"{qualifiers[key]}; {value}" if qualifiers[key] else value
    else:
        qualifiers.setdefault(key, value)


def read_genbank(path: Path) -> list[GenBankRecord]:
    """The records of a GenBank flat file (gzipped or not): names, features and sequences."""
    try:
        with open_text(path) as fh:
            return list(_genbank_records(fh))
    except DECOMPRESSION_ERRORS as exc:
        raise BaconError(f"{path}: truncated or corrupt compressed file ({exc})") from None


@dataclass
class _Pending:
    """A feature being read: its key, location text and qualifiers so far."""

    type: str
    location: str
    qualifiers: dict[str, str] = field(default_factory=dict)
    current: str = ""  # The qualifier being read
    quoted: bool = False  # Inside a quoted value that continues on the next line


def _genbank_records(fh) -> list[GenBankRecord]:  # noqa: C901 - a line-oriented state machine
    records: list[GenBankRecord] = []
    rec: GenBankRecord | None = None
    section = ""
    keyword = ""
    definition: list[str] = []
    seq_chunks: list[str] = []
    pending: _Pending | None = None

    def close_feature() -> None:
        nonlocal pending
        if pending is None or rec is None:
            return
        try:
            location = parse_location(pending.location)
        except ValueError:
            rec.skipped += 1
        else:
            if location is not None:
                rec.features.append(RawFeature(pending.type, "", location, pending.qualifiers))
        pending = None

    for raw in fh:
        line = raw.rstrip("\r\n")
        if line.startswith("LOCUS"):
            fields = line.split()
            length = int(fields[2]) if len(fields) > 3 and fields[2].isdigit() and fields[3] == "bp" else 0
            rec = GenBankRecord(fields[1] if len(fields) > 1 else "", "", length, "circular" in fields, [], "")
            section, keyword, definition, seq_chunks = "header", "", [], []
            continue
        if rec is None:
            continue
        if line.startswith("//"):
            close_feature()
            rec.definition = " ".join(" ".join(definition).split()).removesuffix(".")
            rec.seq = "".join(seq_chunks)
            if rec.length == 0:
                rec.length = len(rec.seq)
            for f in rec.features:
                f.seq = rec.name
            records.append(rec)
            rec = None
            continue
        if section == "sequence":
            seq_chunks.append("".join(ch for ch in line if ch.isalpha()))
            continue
        if line[:1] not in (" ", ""):  # A keyword line
            close_feature()
            keyword, value = line[:12].strip(), line[12:].strip()
            section = "features" if keyword == "FEATURES" else "sequence" if keyword == "ORIGIN" else "header"
            if keyword == "DEFINITION":
                definition = [value]
            elif keyword == "VERSION" and value:
                rec.name = value.split()[0]
            continue
        if section == "header":
            if keyword == "DEFINITION":
                definition.append(line.strip())
            continue
        if section != "features":
            continue
        key = line[5:21].strip()
        if line[:5] == "     " and key and not line.startswith(" " * 21):  # A new feature
            close_feature()
            pending = _Pending(key, line[21:].strip())
            continue
        if pending is None:
            continue
        text = line.strip()
        q = pending.qualifiers
        if pending.quoted:  # A quoted value continues
            value = text.removesuffix('"')
            pending.quoted = not text.endswith('"')
            q[pending.current] = (q[pending.current] + value if pending.current == "translation"
                                  else f"{q[pending.current]} {value}".strip())
        elif text.startswith("/"):
            key, _, value = text[1:].partition("=")
            if value.startswith('"'):
                if len(value) > 1 and value.endswith('"'):
                    value = value[1:-1]
                else:
                    value = value[1:]
                    pending.quoted = True
            _add_qualifier(q, key, value)
            pending.current = key.lower()
        elif not pending.current:  # The location continues
            pending.location += text
        else:  # An unquoted value continues
            q[pending.current] += text
    close_feature()
    return records


def genbank_fasta_records(path: Path) -> list[Record]:
    """The sequences of a GenBank file, named as NCBI's fasta of the same records (for `-r file.gb`)."""
    records = read_genbank(path)
    if not records:
        raise BaconError(f"{path}: no GenBank record (no LOCUS line)")
    missing = [r.name for r in records if not r.seq]
    if missing:
        raise BaconError(f"{path}: no sequence (ORIGIN) for {', '.join(missing[:3])}: give the reference as fasta "
                         "and this file with --annotation")
    return [Record(r.header, r.seq) for r in records]


def read_gff3(path: Path) -> tuple[list[RawFeature], dict[str, int], set[str]]:
    """The features of a GFF3 file (gzipped or not; a ##FASTA section is ignored), the sequence lengths given by
    ##sequence-region lines, and the names of the sequences flagged circular. Multi-line features (the parts of
    a CDS, sharing an ID) are merged. Attribute keys are lower-cased."""
    features: dict[str, RawFeature] = {}  # By ID
    order: list[RawFeature] = []
    lengths: dict[str, int] = {}
    circular: set[str] = set()
    try:
        with open_text(path) as fh:
            for line in fh:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("##sequence-region"):
                    fields = line.split()
                    if len(fields) >= 4 and fields[3].isdigit():
                        lengths[fields[1]] = int(fields[3])
                    continue
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\r\n").split("\t")
                if len(cols) < 8 or not (cols[3].isdigit() and cols[4].isdigit()):
                    continue
                start, end = int(cols[3]), int(cols[4])
                if end < start or start < 1:
                    continue
                attributes: dict[str, str] = {}
                for item in (cols[8] if len(cols) > 8 else "").split(";"):
                    key, eq, value = item.strip().partition("=")
                    if key:
                        _add_qualifier(attributes, unquote(key), unquote(value) if eq else "")
                if attributes.get("is_circular", "").lower() == "true":
                    circular.add(cols[0])
                strand = -1 if cols[6] == "-" else 1
                phase = int(cols[7]) if cols[7] in ("0", "1", "2") else None
                ident = attributes.get("id", "")
                partial_low = "start_range" in attributes
                partial_high = "end_range" in attributes
                same = features.get(ident) if ident else None
                if same is not None and same.type == cols[2] and same.seq == cols[0]:  # Another part
                    same.location.parts = sorted(same.location.parts + [(start, end)])
                    same.location.partial_low |= partial_low
                    same.location.partial_high |= partial_high
                    if phase is not None:
                        same.phases[start] = phase
                    continue
                feature = RawFeature(cols[2], cols[0], Location(strand, [(start, end)], partial_low, partial_high),
                                     attributes, ident,
                                     tuple(p for p in attributes.get("parent", "").split(",") if p))
                if phase is not None:
                    feature.phases[start] = phase
                if ident and ident not in features:
                    features[ident] = feature
                order.append(feature)
    except DECOMPRESSION_ERRORS as exc:
        raise BaconError(f"{path}: truncated or corrupt compressed file ({exc})") from None
    return order, lengths, circular


# ---------------------------------------------------------------------------------------------------------------
# Genes and regions
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class Cds:
    strand: int
    parts: list[tuple[int, int]]  # In the order of translation (5' to 3' on the coding strand)
    codon_start: int = 1
    transl_table: int = DEFAULT_TABLE
    partial5: bool = False  # The start of the coding sequence is missing
    partial3: bool = False
    trans_spliced: bool = False
    strands: list[int] | None = None  # The strand of each part, when they differ (trans-splicing across strands)
    _coding: tuple[int, str] | None = field(default=None, repr=False, compare=False)  # Cache: (id(seq), coding)

    def part_strand(self, k: int) -> int:
        return self.strands[k] if self.strands else self.strand

    def coding_sequence(self, seq: str) -> str:
        """The coding sequence read from `seq` (the reference, 5' to 3' on the coding strand), from the first base
        given by codon_start."""
        if self._coding is None or self._coding[0] != id(seq):
            text = "".join(reverse_complement(seq[s - 1:e]) if self.part_strand(k) < 0 else seq[s - 1:e]
                           for k, (s, e) in enumerate(self.parts))
            self._coding = (id(seq), text[self.codon_start - 1:])
        return self._coding[1]

    def index(self, pos: int) -> tuple[int, int] | None:
        """0-based position of a reference base in the coding sequence (after codon_start) and the strand of the
        part it is in; None outside the coding sequence."""
        offset = 0
        for k, (s, e) in enumerate(self.parts):
            strand = self.part_strand(k)
            if s <= pos <= e:
                return offset + (e - pos if strand < 0 else pos - s) - (self.codon_start - 1), strand
            offset += e - s + 1
        return None


@dataclass
class Gene:
    seq: str
    name: str
    kind: str  # "CDS", "tRNA", "rRNA", "ncRNA", "pseudogene" or "other"
    strand: int
    extent: list[tuple[int, int]]  # The gene feature's ranges (the hull of the parts without a gene feature)
    cds: list[Cds] = field(default_factory=list)
    exons: list[tuple[int, int]] = field(default_factory=list)  # Parts of the CDS, tRNA and rRNA features
    product: str = ""
    pseudo: bool = False
    snps: int = 0  # SNP positions inside the gene, counted by annotate_snps

    @property
    def start(self) -> int:
        return self.extent[0][0]

    @property
    def end(self) -> int:
        return self.extent[-1][1]

    @property
    def length(self) -> int:
        return sum(e - s + 1 for s, e in self.extent)

    def contains(self, pos: int) -> bool:
        return any(s <= pos <= e for s, e in self.extent)

    def context(self, pos: int) -> str:
        """What the position is in: CDS, tRNA, rRNA, ncRNA, pseudogene, intron, UTR or gene."""
        if self.kind == "pseudogene":
            return "pseudogene"
        if any(s <= pos <= e for s, e in self.exons):
            return self.kind if self.kind != "other" else "gene"
        if self.exons and self.exons[0][0] < pos < self.exons[-1][1]:
            return "intron"
        return "UTR" if self.exons else "gene"

    def kind_text(self) -> str:
        return {"CDS": "protein-coding gene", "tRNA": "tRNA gene", "rRNA": "rRNA gene", "ncRNA": "non-coding RNA",
                "pseudogene": "pseudogene"}.get(self.kind, "gene")


@dataclass
class Region:
    name: str
    start: int
    end: int  # Smaller than start when the region spans the origin of a circular sequence
    length: int

    def contains(self, pos: int) -> bool:
        if self.start <= self.end:
            return self.start <= pos <= self.end
        return pos >= self.start or pos <= self.end


def _gene_name(q: dict[str, str], fallback: str) -> str:
    for key in ("gene", "locus_tag", "product", "name", "id"):
        if q.get(key):
            return q[key].split("; ")[0]
    return fallback


def _is_pseudo(q: dict[str, str]) -> bool:
    return ("pseudo" in q and q["pseudo"].lower() != "false") or "pseudogene" in q \
        or q.get("gene_biotype") == "pseudogene"


def _cds(f: RawFeature, default_table: int) -> Cds:
    loc = f.location
    partial5, partial3 = loc.partial_low, loc.partial_high
    if loc.strand < 0:
        partial5, partial3 = partial3, partial5
    codon_start = 1
    if f.qualifiers.get("codon_start", "").strip() in ("1", "2", "3"):
        codon_start = int(f.qualifiers["codon_start"])
    elif f.phases:  # GFF3: the phase of the first part in the direction of translation
        codon_start = f.phases.get(loc.ordered[0][0], 0) + 1
    table = f.qualifiers.get("transl_table", "").strip()
    trans = "trans_splicing" in f.qualifiers or "trans" in f.qualifiers.get("exception", "").lower()
    strands = [strand for _, _, strand in loc.ordered]
    return Cds(loc.strand, [(s, e) for s, e, _ in loc.ordered], codon_start,
               int(table) if table.isdigit() else default_table, partial5, partial3, trans,
               strands if len(set(strands)) > 1 else None)


def _attach(gene: Gene, f: RawFeature, default_table: int) -> None:
    """Add a CDS, tRNA, rRNA or other RNA feature to its gene."""
    if _is_pseudo(f.qualifiers):
        gene.pseudo = True
    if f.type == "CDS":
        gene.cds.append(_cds(f, default_table))
        if gene.kind in ("other", "ncRNA"):
            gene.kind = "CDS"
    elif f.type in RNA_TYPES and gene.kind == "other":
        gene.kind = RNA_TYPES[f.type]
    if f.type == "CDS" or f.type in RNA_TYPES:
        gene.exons = sorted(set(gene.exons) | set(f.location.parts))
    if not gene.product and f.qualifiers.get("product"):
        gene.product = f.qualifiers["product"].split("; ")[0]


def _new_gene(f: RawFeature, default_table: int) -> Gene:
    """A gene made from a CDS or RNA feature without a gene feature of its own."""
    loc = f.location
    trans = "trans_splicing" in f.qualifiers or "trans" in f.qualifiers.get("exception", "").lower()
    extent = list(loc.parts) if trans else [(loc.start, loc.end)]
    gene = Gene(f.seq, _gene_name(f.qualifiers, f"{f.type}:{loc.start}"), "other", loc.strand, extent,
                pseudo=_is_pseudo(f.qualifiers))
    _attach(gene, f, default_table)
    return gene


def _finish(genes: list[Gene]) -> list[Gene]:
    for g in genes:
        if g.pseudo:
            g.kind = "pseudogene"
        g.exons.sort()
    return sorted(genes, key=lambda g: (g.start, g.end))


def _overlaps(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> bool:
    return any(s1 <= e2 and s2 <= e1 for s1, e1 in a for s2, e2 in b)


def genes_from_genbank(features: list[RawFeature], default_table: int = DEFAULT_TABLE) -> list[Gene]:
    """Group GenBank features into genes: a CDS, tRNA or rRNA joins the gene feature with the same /gene or
    /locus_tag that overlaps it (or, without a name, any overlapping gene feature on its strand); without one, it
    is a gene of its own."""
    genes: list[Gene] = []
    by_key: dict[tuple[str, str], list[Gene]] = {}
    for f in features:
        if f.type == "gene":
            gene = Gene(f.seq, _gene_name(f.qualifiers, f"gene:{f.location.start}"), "other", f.location.strand,
                        list(f.location.parts), pseudo=_is_pseudo(f.qualifiers))
            genes.append(gene)
            for key in ("gene", "locus_tag"):
                if f.qualifiers.get(key):
                    by_key.setdefault((f.seq, f.qualifiers[key]), []).append(gene)
    for f in features:
        if f.type != "CDS" and f.type not in RNA_TYPES:
            continue
        keys = [(f.seq, f.qualifiers[k]) for k in ("locus_tag", "gene") if f.qualifiers.get(k)]
        candidates = [g for key in keys for g in by_key.get(key, []) if _overlaps(g.extent, f.location.parts)]
        if not candidates and not keys:
            candidates = [g for g in genes if g.seq == f.seq and g.strand == f.location.strand
                          and _overlaps(g.extent, f.location.parts)]
        if candidates:
            _attach(candidates[0], f, default_table)
        else:
            genes.append(_new_gene(f, default_table))
    return _finish(genes)


def genes_from_gff3(features: list[RawFeature], default_table: int = DEFAULT_TABLE) -> list[Gene]:
    """Group GFF3 features into genes through their Parent links (gene > mRNA/tRNA/rRNA > CDS/exon)."""
    by_id = {f.id: f for f in features if f.id}
    gene_of: dict[str, Gene] = {}
    genes: list[Gene] = []
    for f in features:
        if f.type in ("gene", "pseudogene"):
            gene = Gene(f.seq, _gene_name(f.qualifiers, f"gene:{f.location.start}"), "other", f.location.strand,
                        list(f.location.parts), pseudo=f.type == "pseudogene" or _is_pseudo(f.qualifiers))
            genes.append(gene)
            if f.id:
                gene_of[f.id] = gene

    def ancestor(f: RawFeature) -> Gene | None:
        seen, queue = set(), list(f.parents)
        while queue:
            parent = queue.pop(0)
            if parent in gene_of:
                return gene_of[parent]
            if parent in by_id and parent not in seen:
                seen.add(parent)
                queue += by_id[parent].parents
        return None

    for f in features:
        if f.type == "CDS" or f.type in RNA_TYPES:
            gene = ancestor(f)
            if gene is None:
                genes.append(_new_gene(f, default_table))
            else:
                _attach(gene, f, default_table)
        elif f.type == "exon":
            gene = ancestor(f)
            if gene is not None and gene.kind != "CDS":
                gene.exons = sorted(set(gene.exons) | set(f.location.parts))
    return _finish(genes)


_REGION = re.compile(r"^(?:(?P<lsc>LSC|large single[ -]copy)|(?P<ssc>SSC|small single[ -]copy)"
                     r"|(?P<ir>IR|inverted[ -]repeat)(?:[ _-]?(?P<copy>[AB]))?)(?: region)?$", re.I)


_IR_COPY = re.compile(r"\b(?:IR|inverted[ -]repeats?(?: region)?)[ _-]?(?:IR)?([AB])\b", re.I)


def region_label(f: RawFeature) -> str | None:
    """'LSC', 'SSC', 'IRa', 'IRb' or 'IR' for a feature annotating a region of a plastome, else None."""
    if f.type not in REGION_TYPES:
        return None
    texts: list[str] = []
    for key in ("note", "standard_name", "name", "rpt_family", "gene", "product", "rpt_type"):
        if f.qualifiers.get(key):
            texts += re.split(r"[;,]", f.qualifiers[key])
    for text in texts:  # A copy named anywhere: "IRa", "inverted repeat B", "inverted repeat region IRb"
        copy = _IR_COPY.search(text)
        if copy:
            return "IR" + copy.group(1).lower()
    for text in texts:
        match = _REGION.match(text.strip())
        if match:
            if match.group("lsc"):
                return "LSC"
            if match.group("ssc"):
                return "SSC"
            return "IR" + (match.group("copy") or "").lower()
    if f.type == "repeat_region" and f.qualifiers.get("rpt_type", "").lower() == "inverted":
        return "IR"
    return None


def derive_regions(features: list[RawFeature], length: int) -> list[Region]:
    """The LSC/IRb/SSC/IRa band of a plastome from its annotated inverted repeats: the two repeats (named as the
    annotation names them, else by convention: IRb follows the LSC), and the single-copy regions as the gaps
    between them, the larger one being the LSC. One of them may span the origin. Without two inverted repeats
    (or with repeats that overlap or touch), no regions."""
    irs: list[tuple[int, int, str]] = []
    for f in features:
        label = region_label(f)
        if label and label.startswith("IR"):
            irs += [(s, e, label) for s, e in f.location.parts if e - s + 1 >= MIN_INVERTED_REPEAT and e <= length]
    if len(irs) != 2:
        return []
    irs.sort()
    (a_start, a_end, a_label), (b_start, b_end, b_label) = irs
    gap_a = (a_end + 1, b_start - 1)  # Between the repeats
    gap_b = (b_end % length + 1, a_start - 1 if a_start > 1 else length)  # Around the origin
    len_a = gap_a[1] - gap_a[0] + 1
    len_b = (length - b_end) + (a_start - 1)
    if len_a <= 0 or len_b <= 0:
        return []
    if {a_label, b_label} != {"IRa", "IRb"}:
        a_label, b_label = ("IRa", "IRb") if len_a >= len_b else ("IRb", "IRa")
    lsc_first = len_a >= len_b
    regions = [Region(a_label, a_start, a_end, a_end - a_start + 1),
               Region("LSC" if lsc_first else "SSC", *gap_a, len_a),
               Region(b_label, b_start, b_end, b_end - b_start + 1),
               Region("SSC" if lsc_first else "LSC", *gap_b, len_b)]
    return sorted(regions, key=lambda r: r.start)


# ---------------------------------------------------------------------------------------------------------------
# The annotation of a reference
# ---------------------------------------------------------------------------------------------------------------

class SequenceAnnotation:
    """The genes and regions of one reference sequence, indexed by position."""

    def __init__(self, name: str, length: int, genes: list[Gene], regions: list[Region], circular: bool):
        self.name, self.length, self.circular = name, length, circular
        self.genes = sorted(genes, key=lambda g: (g.start, g.end))
        self.regions = regions
        self._starts = [g.start for g in self.genes]
        self._by_end = sorted(self.genes, key=lambda g: (g.end, g.start))
        self._ends = [g.end for g in self._by_end]
        self._bins: dict[int, list[Gene]] = {}
        for g in self.genes:
            for s, e in g.extent:
                for b in range(s // _BIN, e // _BIN + 1):
                    self._bins.setdefault(b, []).append(g)

    def genes_at(self, pos: int) -> list[Gene]:
        seen: set[int] = set()
        found = []
        for g in self._bins.get(pos // _BIN, []):
            if id(g) not in seen and g.contains(pos):
                seen.add(id(g))
                found.append(g)
        return found

    def neighbours(self, pos: int) -> tuple[Gene | None, Gene | None]:
        """The nearest gene ending before `pos` and the nearest starting after it (around the origin when the
        sequence is circular)."""
        if not self.genes:
            return None, None
        i = bisect_left(self._ends, pos)
        before = self._by_end[i - 1] if i > 0 else (self._by_end[-1] if self.circular else None)
        j = bisect_right(self._starts, pos)
        after = self.genes[j] if j < len(self.genes) else (self.genes[0] if self.circular else None)
        return before, after

    def region_at(self, pos: int) -> str:
        for r in self.regions:
            if r.contains(pos):
                return r.name
        return ""


@dataclass
class Annotation:
    file: Path
    format: str
    sequences: dict[str, SequenceAnnotation]  # By reference sequence name, those with features
    warnings: list[str] = field(default_factory=list)
    skipped: int = 0  # Features ignored (unparsable location, outside the sequence)

    @property
    def genes(self) -> int:
        return sum(len(s.genes) for s in self.sequences.values())

    @property
    def tables(self) -> list[int]:
        return sorted({c.transl_table for s in self.sequences.values() for g in s.genes for c in g.cds})

    @property
    def has_regions(self) -> bool:
        return any(s.regions for s in self.sequences.values())


def load_annotation(path: Path, sequences: list[tuple[str, int]],
                    default_table: int = DEFAULT_TABLE) -> Annotation:
    """Read a GenBank or GFF3 annotation and match it to the reference sequences (name, length). An unreadable
    file or an unknown format is a BaconError; annotated sequences that match no reference sequence, and
    features beyond the end of their sequence, are left out with a warning (in `warnings`). A single annotated
    sequence is taken for a single reference sequence of the same length, whatever its name."""
    if not path.is_file():
        raise BaconError(f"Annotation file not found: {path}")
    fmt = annotation_format(path)
    if fmt is None:
        raise BaconError(f"{path}: not a GenBank (.gb, .gbk) or GFF3 (.gff, .gff3) annotation")
    lengths: dict[str, int] = {}
    circular: set[str] = set()
    skipped = 0
    if fmt == "genbank":
        records = read_genbank(path)
        if not records:
            raise BaconError(f"{path}: no GenBank record (no LOCUS line)")
        features = [f for r in records for f in r.features]
        lengths = {r.name: r.length for r in records}
        circular = {r.name for r in records if r.circular}
        skipped = sum(r.skipped for r in records)
    else:
        features, lengths, circular = read_gff3(path)
        if not features:
            raise BaconError(f"{path}: no feature (not a GFF3 file?)")
    annotated = list(dict.fromkeys([*lengths, *(f.seq for f in features)]))
    reference = dict(sequences)
    rename: dict[str, str] = {name: name for name in annotated if name in reference}
    warnings = []
    if not rename and len(annotated) == 1 and len(sequences) == 1:
        (ref_name, ref_length), = sequences
        if lengths.get(annotated[0], ref_length) == ref_length:
            rename = {annotated[0]: ref_name}
            warnings.append(f"{path.name}: the annotated sequence {annotated[0]!r} is taken for the reference "
                            f"sequence {ref_name!r} (same length, {ref_length:,} bp)")
    unmatched = [name for name in annotated if name not in rename]
    if unmatched:
        text = ", ".join(unmatched[:5]) + (" …" if len(unmatched) > 5 else "")
        warnings.append(f"{path.name}: {len(unmatched)} annotated sequence(s) match no reference sequence by name "
                        f"and are ignored: {text}; the reference's names are "
                        f"{', '.join(name for name, _ in sequences[:5])}")
    kept: dict[str, list[RawFeature]] = {}
    out_of_range = 0
    for f in features:
        name = rename.get(f.seq)
        if name is None:
            continue
        if f.location.end > reference[name]:
            out_of_range += 1
            continue
        f.seq = name
        kept.setdefault(name, []).append(f)
    if out_of_range:
        warnings.append(f"{path.name}: {out_of_range} feature(s) beyond the end of their reference sequence are "
                        "ignored (is it the annotation of this reference?)")
    builder = genes_from_genbank if fmt == "genbank" else genes_from_gff3
    annotation = Annotation(path, fmt, {}, warnings, skipped + out_of_range)
    for name, length in sequences:
        if name in kept:
            is_circular = any(old in circular for old, new in rename.items() if new == name)
            annotation.sequences[name] = SequenceAnnotation(
                name, length, builder(kept[name], default_table), derive_regions(kept[name], length), is_circular)
    if rename and not annotation.sequences:
        warnings.append(f"{path.name}: no feature on the reference sequences")
    return annotation


# ---------------------------------------------------------------------------------------------------------------
# Effects of SNPs
# ---------------------------------------------------------------------------------------------------------------

_CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_BASES = "TCAG"
STANDARD_CODE = {a + b + c: _CODE[16 * i + 4 * j + k]
                 for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}
# Differences from the standard code and the start codons of the known tables (NCBI genetic codes).
TABLES: dict[int, tuple[dict[str, str], set[str]]] = {
    1: ({}, {"ATG"}),
    4: ({"TGA": "W"}, {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA", "TTA"}),
    11: ({}, {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}),
}
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def translate(codon: str, table: int = DEFAULT_TABLE, start: bool = False) -> str:
    """One-letter amino acid of a codon ('*' for stop, 'X' when ambiguous); the initiation codon is M when it is a
    start codon of the table."""
    changes, starts = TABLES.get(table, TABLES[DEFAULT_TABLE])
    codon = codon.upper()
    if start and codon in starts:
        return "M"
    return changes.get(codon) or STANDARD_CODE.get(codon, "X")


@dataclass
class Effect:
    gene: str
    alt: str  # The alternate allele, as in the VCF
    codon: str
    codon_alt: str
    position: int  # 1-based amino-acid position
    aa: str
    aa_alt: str
    kind: str  # synonymous, missense, nonsense, stop lost, stop retained, start lost, start retained

    @property
    def change(self) -> str:
        return f"{self.aa}{self.position}{self.aa_alt}"

    @property
    def codons(self) -> str:
        return f"{self.codon}>{self.codon_alt}"


def cds_effect(cds: Cds, gene: str, pos: int, ref: str, alt: str, seq: str) -> Effect | None:
    """The effect of a SNP on a coding sequence; None when the position is outside it, the codon is incomplete
    (a partial CDS) or contains N, or the VCF's reference base does not match the reference."""
    if len(ref) != 1 or len(alt) != 1 or alt.upper() not in "ACGT":
        return None
    found = cds.index(pos)
    if found is None or found[0] < 0:
        return None
    i, strand = found
    coding = cds.coding_sequence(seq)
    codon_index, within = divmod(i, 3)
    codon = coding[3 * codon_index:3 * codon_index + 3]
    if len(codon) < 3 or "N" in codon:
        return None
    base, alt_base = (ref.upper(), alt.upper())
    if strand < 0:
        base, alt_base = reverse_complement(base), reverse_complement(alt_base)
    if codon[within] != base:
        return None
    codon_alt = codon[:within] + alt_base + codon[within + 1:]
    first = codon_index == 0 and not cds.partial5
    aa = translate(codon, cds.transl_table, start=first)
    aa_alt = translate(codon_alt, cds.transl_table, start=first)
    if first and translate(codon, cds.transl_table, start=True) == "M":
        if aa_alt == "M":
            kind = "start retained"
        else:
            kind, aa_alt = "start lost", translate(codon_alt, cds.transl_table)
    elif aa == aa_alt:
        kind = "stop retained" if aa == "*" else "synonymous"
    elif aa == "*":
        kind = "stop lost"
    elif aa_alt == "*":
        kind = "nonsense"
    else:
        kind = "missense"
    return Effect(gene, alt, codon, codon_alt, codon_index + 1, aa, aa_alt, kind)


@dataclass
class SnpAnnotation:
    region: str  # "" without a region band
    genes: list[Gene]
    context: str  # CDS, intron, tRNA, rRNA, ..., or "intergenic between X and Y"
    effects: list[Effect]


def annotate_snp(seq_ann: SequenceAnnotation, pos: int, ref: str, alts: list[str], seq: str) -> SnpAnnotation:
    genes = seq_ann.genes_at(pos)
    region = seq_ann.region_at(pos)
    effects = []
    if genes:
        contexts = list(dict.fromkeys(g.context(pos) for g in genes))
        for g in genes:
            if g.kind == "pseudogene":
                continue
            for cds in g.cds:
                for alt in alts:
                    effect = cds_effect(cds, g.name, pos, ref, alt, seq)
                    if effect is not None:
                        effects.append(effect)
        return SnpAnnotation(region, genes, " / ".join(contexts), effects)
    before, after = seq_ann.neighbours(pos)
    if before is not None and after is not None:
        context = f"intergenic between {before.name} and {after.name}"
    elif after is not None:
        context = f"intergenic before {after.name}"
    elif before is not None:
        context = f"intergenic after {before.name}"
    else:
        context = "intergenic"
    return SnpAnnotation(region, [], context, [])


def annotate_snps(annotation: Annotation, snps: list[tuple[str, int, str, str]],
                  seqs: dict[str, str]) -> dict[tuple[str, int], SnpAnnotation]:
    """Annotate the SNPs (sequence, position, REF, ALT alleles separated by commas) of the annotated sequences,
    and count the SNP positions of each gene (`Gene.snps`)."""
    result: dict[tuple[str, int], SnpAnnotation] = {}
    for s in annotation.sequences.values():
        for g in s.genes:
            g.snps = 0
    for chrom, pos, ref, alt in snps:
        seq_ann = annotation.sequences.get(chrom)
        if seq_ann is None or (chrom, pos) in result or chrom not in seqs:
            continue
        info = annotate_snp(seq_ann, pos, ref, [a for a in alt.split(",") if a], seqs[chrom])
        for g in info.genes:
            g.snps += 1
        result[(chrom, pos)] = info
    return result
