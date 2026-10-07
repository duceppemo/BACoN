"""Reference annotations (GenBank or GFF3): the genes and regions of the reference, and the effect of the VCF's
SNPs on its coding sequences. Used by the report only (standard library only).

Sequence names: a GenBank record is named after its VERSION (accession.version, `NC_008096.2`), or its LOCUS name
without one, as NCBI's fasta of the same record; a GFF3 feature after its first column. Translation: each CDS's
`/transl_table` (GFF3: `transl_table=`), else table 11 (bacterial and plastid code); the genetic codes 1, 2, 3, 4,
5, 9, 11, 13 and 14 are known, other tables are translated with the standard code and table 11's start codons.
"""

from __future__ import annotations

import logging
import re
from bisect import bisect_left, bisect_right
from collections import Counter
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
MAX_MAP_PRODUCT = 12  # A gene without a symbol is labelled on the map by its product when it is this short
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


def _join_wrapped(previous: str, text: str) -> str:
    """A quoted value continued on the next line, joined with a space, except where the flat file broke a word:
    after a hyphen (`ATP-` / `dependent`) or a comma followed by a digit (`2,` / `6-diaminopimelate`)."""
    if (previous.endswith("-") and text[:1].isalnum()) or (previous.endswith(",") and text[:1].isdigit()):
        return previous + text
    return f"{previous} {text}".strip()


def _closes_quote(text: str) -> bool:
    """Whether a line ends a quoted GenBank value: it ends with an odd number of quotes (`""` is a quote)."""
    return (len(text) - len(text.rstrip('"'))) % 2 == 1


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
            closed = _closes_quote(text)
            value = (text[:-1] if closed else text).replace('""', '"')
            pending.quoted = not closed
            q[pending.current] = (q[pending.current] + value if pending.current == "translation"
                                  else _join_wrapped(q[pending.current], value))
        elif text.startswith("/"):
            key, _, value = text[1:].partition("=")
            if value.startswith('"'):
                value = value[1:]
                closed = _closes_quote(value)
                value = (value[:-1] if closed else value).replace('""', '"')
                pending.quoted = not closed
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


def _order_parts(loc: Location, numbers: list[int | None], length: int) -> None:
    """Put the parts of a multi-line GFF3 feature in the order of translation: by their `part=` numbers (NCBI)
    when every line has one, else as listed (NCBI lists them 5' to 3', which puts a trans-spliced CDS, or one
    across the origin, out of coordinate order), except for parts listed by ascending coordinate (Ensembl, a
    sorted file): those of a reverse-strand feature are read in descending order, and those of a feature across
    the origin (spanning more than half the sequence, from its first base to its last) start at the high part on
    the + strand and at the low parts, read downwards, on the - strand. Without the sequence's length (no
    ##sequence-region, and no reference sequence of that name) a feature across the origin cannot be told apart.
    The strand of the feature is that of most of its bases."""
    listed = loc.listed or []
    strands = {strand for _, _, strand in listed}
    if all(n is not None for n in numbers) and len(set(numbers)) == len(numbers):
        listed = [part for _, part in sorted(zip(numbers, listed))]
    elif len(strands) == 1 and listed == sorted(listed):
        strand = strands.pop()
        across = bool(length) and listed[-1][1] - listed[0][0] + 1 > length / 2
        if across and listed[0][0] == 1 and listed[-1][1] == length:
            # The parts either side of the origin are either side of the widest gap
            k = max(range(1, len(listed)), key=lambda i: listed[i][0] - listed[i - 1][1])
            low, high = listed[:k], listed[k:]
            listed = high + low if strand > 0 else low[::-1] + high[::-1]
        elif strand < 0 and not across:
            listed = listed[::-1]
    loc.listed = listed
    loc.strand = 1 if sum((e - s + 1) * strand for s, e, strand in listed) >= 0 else -1


def read_gff3(path: Path, known_lengths: dict[str, int] | None = None
              ) -> tuple[list[RawFeature], dict[str, int], set[str]]:
    """The features of a GFF3 file (gzipped or not; a ##FASTA section is ignored), the sequence lengths given by
    ##sequence-region lines, and the names of the sequences flagged circular (`Is_circular=true`). Multi-line
    features (the parts of a CDS, sharing an ID, or CDS lines without an ID sharing a Parent) are merged, their
    parts kept in the order of translation (`Location.listed`); a line of a circular sequence ending beyond its
    length (Bakta writes a feature across the origin so, with end = its end + the length) is split into its two
    parts. The lengths come from the file, else from `known_lengths` (the reference's). Attribute keys are
    lower-cased."""
    features: dict[str, RawFeature] = {}  # By ID (a CDS without one: by sequence and Parent)
    part_numbers: dict[str, list[int | None]] = {}  # The `part=` attribute of each line of a merged feature
    order: list[RawFeature] = []
    lengths: dict[str, int] = {}
    circular: set[str] = set()
    known_lengths = known_lengths or {}
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
                parents = tuple(p for p in attributes.get("parent", "").split(",") if p)
                key = ident or (f"\t{cols[0]}\t{attributes['parent']}" if cols[2] == "CDS" and parents else "")
                partial_low = "start_range" in attributes
                partial_high = "end_range" in attributes
                part = attributes.get("part", "")
                number = int(part) if part.isdigit() else None
                length = lengths.get(cols[0]) or known_lengths.get(cols[0], 0)
                pieces = [(start, end)]  # In the order of translation
                if length and end > length >= start and end - length < start and cols[0] in circular:
                    pieces = [(start, length), (1, end - length)]  # Across the origin, as Bakta writes it
                    if strand < 0:
                        pieces.reverse()
                same = features.get(key) if key else None
                if same is not None and same.type == cols[2] and same.seq == cols[0]:  # Another part
                    same.location.listed = [*(same.location.listed or []), *((s, e, strand) for s, e in pieces)]
                    same.location.parts = sorted(same.location.parts + pieces)
                    same.location.partial_low |= partial_low
                    same.location.partial_high |= partial_high
                    if phase is not None:
                        same.phases[pieces[0][0]] = phase
                    part_numbers[key].append(number)
                    continue
                location = Location(strand, sorted(pieces), partial_low, partial_high,
                                    [(s, e, strand) for s, e in pieces])
                feature = RawFeature(cols[2], cols[0], location, attributes, ident, parents)
                if phase is not None:
                    feature.phases[pieces[0][0]] = phase
                if key and key not in features:
                    features[key] = feature
                    part_numbers[key] = [number]
                order.append(feature)
    except DECOMPRESSION_ERRORS as exc:
        raise BaconError(f"{path}: truncated or corrupt compressed file ({exc})") from None
    for key, numbers in part_numbers.items():
        if len(numbers) > 1:
            seq = features[key].seq
            _order_parts(features[key].location, numbers, lengths.get(seq) or known_lengths.get(seq, 0))
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
    transl_except: list[tuple[int, int]] = field(default_factory=list)  # Codons read otherwise (/transl_except)
    _coding: tuple[int, str] | None = field(default=None, repr=False, compare=False)  # Cache: (id(seq), coding)

    def part_strand(self, k: int) -> int:
        return self.strands[k] if self.strands else self.strand

    def excepted(self, pos: int) -> bool:
        """Whether the position is in a codon with a translational exception (selenocysteine, pyrrolysine, an
        edited or a partial stop codon)."""
        return any(s <= pos <= e for s, e in self.transl_except)

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
    symbol: str = ""  # The gene symbol (/gene; GFF3 gene=, or a Name that is not an identifier); "" without one

    @property
    def start(self) -> int:
        return self.extent[0][0]

    @property
    def label(self) -> str:
        """The gene as the report shows it: its name when that is its gene symbol, else the identifier followed
        by the product (or the symbol) in parentheses: `LK299_pgr007 (23S ribosomal RNA)`."""
        detail = self.symbol or self.product
        if self.symbol == self.name or not detail or detail == self.name:
            return self.name
        return f"{self.name} ({detail})"

    @property
    def map_label(self) -> str:
        """A short label for the genome map: the name, or, for a gene without a symbol, its product when it is
        short (at most MAX_MAP_PRODUCT characters: `tRNA-Val` rather than `OrsajCt141`)."""
        if self.symbol or not self.product or len(self.product) > MAX_MAP_PRODUCT:
            return self.name
        return self.product

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
    for key in ("gene", "locus_tag", "product", "name", "gene_id", "id"):  # gene_id: Ensembl, before `gene:` IDs
        if q.get(key):
            return q[key].split("; ")[0]
    return fallback


_LOCUS_TAG = re.compile(r"^[A-Za-z][A-Za-z0-9]*_[A-Za-z0-9]+$")  # PREFIX_number, as NCBI assigns them


def _symbol(q: dict[str, str], gff3_gene: bool = False) -> str:
    """The gene symbol of a feature: its `gene` qualifier (GenBank /gene, GFF3 gene=); for a GFF3 gene feature,
    also its `Name` unless that is an identifier: the locus tag, the product, the gene_id or the ID, or the ID
    without its type prefix (`gene-`, `gene:`) when the feature has a locus tag or the Name looks like one (NCBI
    names a gene without a symbol after its locus tag; a third party's `ID=gene-matK;Name=matK` is a symbol).
    Nothing else: a locus tag or a product is not a symbol, and no symbol is made up from a product."""
    if q.get("gene"):
        return q["gene"].split("; ")[0]
    if not gff3_gene:
        return ""
    name = q.get("name", "").split("; ")[0]
    ident = q.get("id", "")
    if not name or name in (q.get("locus_tag"), q.get("product"), q.get("gene_id"), ident):
        return ""
    if ident.endswith(("-" + name, ":" + name)) and ("locus_tag" in q or _LOCUS_TAG.match(name)):
        return ""
    return name


def _is_pseudo(q: dict[str, str]) -> bool:
    return ("pseudo" in q and q["pseudo"].lower() != "false") or "pseudogene" in q \
        or any(q.get(key, "").lower().endswith("pseudogene") for key in ("gene_biotype", "biotype"))


_TRANSL_EXCEPT = re.compile(r"pos:\s*((?:complement\()?[<>]?\d+(?:\.\.[<>]?\d+)?\)?)\s*,\s*aa:\s*\w+", re.I)


def _transl_except(text: str) -> list[tuple[int, int]]:
    """The reference ranges of `/transl_except=(pos:a..b,aa:Sec)` qualifiers (several joined by '; ' or ',')."""
    ranges: list[tuple[int, int]] = []
    for match in _TRANSL_EXCEPT.finditer(text):
        try:
            location = parse_location(match.group(1))
        except ValueError:
            continue
        if location is not None:
            ranges += location.parts
    return ranges


def _is_trans_spliced(f: RawFeature) -> bool:
    return "trans_splicing" in f.qualifiers or "trans" in f.qualifiers.get("exception", "").lower()


def _spliced_apart(loc: Location, length: int) -> bool:
    """Whether the parts of a location are not one stretch of the sequence read in order: parts on both strands,
    a part followed by one upstream of it (trans-splicing, or a feature across the origin listed 5' to 3'), or
    parts spanning more than half the sequence (across the origin)."""
    ordered = loc.ordered
    if len(ordered) < 2:
        return False
    if len({strand for _, _, strand in ordered}) > 1:
        return True
    for (s1, e1, strand), (s2, e2, _) in zip(ordered, ordered[1:]):
        if (strand > 0 and s2 <= e1) or (strand < 0 and e2 >= s1):
            return True
    return bool(length) and loc.end - loc.start + 1 > length / 2


def _cds(f: RawFeature, default_table: int, length: int = 0) -> Cds:
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
    trans = _is_trans_spliced(f) or _spliced_apart(loc, length)
    strands = [strand for _, _, strand in loc.ordered]
    return Cds(loc.strand, [(s, e) for s, e, _ in loc.ordered], codon_start,
               int(table) if table.isdigit() else default_table, partial5, partial3, trans,
               strands if len(set(strands)) > 1 else None, _transl_except(f.qualifiers.get("transl_except", "")))


def _attach(gene: Gene, f: RawFeature, default_table: int, length: int = 0) -> None:
    """Add a CDS, tRNA, rRNA or other RNA feature to its gene."""
    if _is_pseudo(f.qualifiers):
        gene.pseudo = True
    if f.type == "CDS":
        gene.cds.append(_cds(f, default_table, length))
        if gene.kind in ("other", "ncRNA"):
            gene.kind = "CDS"
    elif f.type in RNA_TYPES and gene.kind == "other":
        gene.kind = RNA_TYPES[f.type]
    if f.type == "CDS" or f.type in RNA_TYPES:
        gene.exons = sorted(set(gene.exons) | set(f.location.parts))
    if not gene.product and f.qualifiers.get("product"):
        gene.product = f.qualifiers["product"].split("; ")[0]
    if not gene.symbol:
        gene.symbol = _symbol(f.qualifiers)


def _new_gene(f: RawFeature, default_table: int, length: int = 0) -> Gene:
    """A gene made from a CDS or RNA feature without a gene feature of its own: its extent is the hull of the
    parts, or the parts themselves when they are not one stretch of the sequence (trans-splicing, a feature
    across the origin)."""
    loc = f.location
    apart = _is_trans_spliced(f) or _spliced_apart(loc, length)
    extent = list(loc.parts) if apart else [(loc.start, loc.end)]
    gene = Gene(f.seq, _gene_name(f.qualifiers, f"{f.type}:{loc.start}"), "other", loc.strand, extent,
                pseudo=_is_pseudo(f.qualifiers))
    _attach(gene, f, default_table, length)
    return gene


def _finish(genes: list[Gene]) -> list[Gene]:
    for g in genes:
        if g.pseudo:
            g.kind = "pseudogene"
        g.exons.sort()
    return sorted(genes, key=lambda g: (g.start, g.end))


def _overlaps(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> bool:
    return any(s1 <= e2 and s2 <= e1 for s1, e1 in a for s2, e2 in b)


def genes_from_genbank(features: list[RawFeature], default_table: int = DEFAULT_TABLE,
                       length: int = 0) -> list[Gene]:
    """Group GenBank features into genes: a CDS, tRNA or rRNA joins the gene feature with the same /gene or
    /locus_tag that overlaps it (or, without a name, any overlapping gene feature on its strand); without one, it
    is a gene of its own. `length`, the sequence's, tells a feature across the origin from a spliced one."""
    genes: list[Gene] = []
    by_key: dict[tuple[str, str], list[Gene]] = {}
    for f in features:
        if f.type == "gene":
            gene = Gene(f.seq, _gene_name(f.qualifiers, f"gene:{f.location.start}"), "other", f.location.strand,
                        list(f.location.parts), pseudo=_is_pseudo(f.qualifiers), symbol=_symbol(f.qualifiers))
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
            _attach(candidates[0], f, default_table, length)
        else:
            genes.append(_new_gene(f, default_table, length))
    return _finish(genes)


GFF3_GENE_TYPES = ("gene", "pseudogene", "ncRNA_gene")  # ncRNA_gene: Ensembl's tRNA and rRNA genes


def genes_from_gff3(features: list[RawFeature], default_table: int = DEFAULT_TABLE, length: int = 0) -> list[Gene]:
    """Group GFF3 features into genes through their Parent links (gene > mRNA/tRNA/rRNA > CDS/exon). A CDS or
    RNA feature without a Parent joins the gene feature on its strand that contains it and has its locus_tag
    or gene (or, when it has neither, any gene feature containing it: NCBI writes the rRNA of a gene with
    `gene_biotype=other` without a Parent), else it is a gene of its own."""
    by_id = {f.id: f for f in features if f.id}
    gene_of: dict[str, Gene] = {}
    genes: list[Gene] = []
    gene_features: list[tuple[RawFeature, Gene]] = []
    for f in features:
        if f.type in GFF3_GENE_TYPES:
            gene = Gene(f.seq, _gene_name(f.qualifiers, f"gene:{f.location.start}"), "other", f.location.strand,
                        list(f.location.parts), pseudo=f.type == "pseudogene" or _is_pseudo(f.qualifiers),
                        symbol=_symbol(f.qualifiers, gff3_gene=True))
            genes.append(gene)
            gene_features.append((f, gene))
            if f.id:
                gene_of[f.id] = gene

    def enclosing(f: RawFeature) -> Gene | None:
        keys = {f.qualifiers[k] for k in ("locus_tag", "gene") if f.qualifiers.get(k)}
        for raw, gene in gene_features:
            if raw.seq != f.seq or raw.location.strand != f.location.strand \
                    or not all(any(s <= ps and pe <= e for s, e in gene.extent) for ps, pe in f.location.parts):
                continue
            own = {raw.qualifiers[k] for k in ("locus_tag", "gene") if raw.qualifiers.get(k)}
            if not keys or keys & own:
                return gene
        return None

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
            gene = ancestor(f) or (enclosing(f) if not f.parents else None)
            if gene is None:
                genes.append(_new_gene(f, default_table, length))
            else:
                _attach(gene, f, default_table, length)
        elif f.type == "exon":
            gene = ancestor(f)
            if gene is not None and gene.kind != "CDS":
                gene.exons = sorted(set(gene.exons) | set(f.location.parts))
    return _finish(genes)


_REGION = re.compile(r"^(?:(?P<lsc>LSC|large single[ -]copy)|(?P<ssc>SSC|small single[ -]copy)"
                     r"|(?P<ir>IR|inverted[ -]repeat)(?:[ _-]?(?P<copy>[AB]))?)(?: region)?$", re.I)


# A copy named at the start of a note (or of one of its ';'/','-separated items): "IRa", "inverted repeat B",
# "inverted repeat region IRb"; not a mention further in ("... in IRA", "junction LSC-IRB")
_IR_COPY = re.compile(r"^\s*(?:IR|inverted[ -]repeats?(?: region)?)[ _-]?(?:IR)?([AB])\b", re.I)


def region_label(f: RawFeature) -> str | None:
    """'LSC', 'SSC', 'IRa', 'IRb' or 'IR' for a feature annotating a region of a plastome, else None."""
    if f.type not in REGION_TYPES:
        return None
    texts: list[str] = []
    for key in ("note", "standard_name", "name", "rpt_family", "gene", "product", "rpt_type"):
        if f.qualifiers.get(key):
            texts += re.split(r"[;,]", f.qualifiers[key])
    for text in texts:
        copy = _IR_COPY.match(text)
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


@dataclass
class _Repeat:
    """An annotated inverted repeat: its pieces of the sequence (one, or two when it is given across the origin),
    merged from the features annotating it."""

    pieces: list[tuple[int, int]]
    label: str  # "IRa", "IRb" or "IR"
    inverted: bool = False  # A repeat_region with /rpt_type=inverted

    def touches(self, pieces: list[tuple[int, int]], length: int) -> bool:
        return _overlaps([(max(s - 1, 1), min(e + 1, length)) for s, e in self.pieces], pieces) or (
            any(e == length for _, e in self.pieces) and any(s == 1 for s, _ in pieces)) or (
            any(s == 1 for s, _ in self.pieces) and any(e == length for _, e in pieces))

    def interval(self, length: int) -> tuple[int, int, int]:
        """(start, end, length), the end smaller than the start when the repeat spans the origin."""
        merged: list[tuple[int, int]] = []
        for s, e in sorted(self.pieces):
            if merged and s <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        if len(merged) == 2 and merged[0][0] == 1 and merged[1][1] == length:
            return merged[1][0], merged[0][1], length - merged[1][0] + 1 + merged[0][1]
        start, end = merged[0][0], merged[-1][1]
        return start, end, end - start + 1


def _inverted_repeats(features: list[RawFeature], length: int) -> list[_Repeat]:
    """The annotated inverted repeats of a sequence, each annotated once: the features labelling the same (or an
    unnamed) repeat over the same stretch are merged; those of less than MIN_INVERTED_REPEAT bp are left out."""
    repeats: list[_Repeat] = []
    for f in features:
        label = region_label(f)
        if not label or not label.startswith("IR"):
            continue
        pieces = [(s, e) for s, e in f.location.parts if e <= length]
        inverted = f.type == "repeat_region" and f.qualifiers.get("rpt_type", "").lower() == "inverted"
        for repeat in repeats:
            if (repeat.label == label or "IR" in (repeat.label, label)) and repeat.touches(pieces, length):
                repeat.pieces += pieces
                repeat.label = label if repeat.label == "IR" else repeat.label
                repeat.inverted |= inverted
                break
        else:
            repeats.append(_Repeat(pieces, label, inverted))
    return [r for r in repeats if r.interval(length)[2] >= MIN_INVERTED_REPEAT]


def derive_regions(features: list[RawFeature], length: int) -> list[Region]:
    """The LSC/IRb/SSC/IRa band of a plastome from its annotated inverted repeats: the two repeats (named as the
    annotation names them, else by convention: IRb follows the LSC), and the single-copy regions as the gaps
    between them, the larger one being the LSC. One region may span the origin (a repeat given as
    join(x..length,1..y) too). With more than two repeats annotated, the two `repeat_region /rpt_type=inverted`,
    else the pair named IRa and IRb, else the two longest; without two (or with repeats that overlap or touch),
    no regions."""
    repeats = _inverted_repeats(features, length)
    if len(repeats) > 2:
        for chosen in ([r for r in repeats if r.inverted], [r for r in repeats if r.label != "IR"],
                       sorted(repeats, key=lambda r: -r.interval(length)[2])[:2]):
            if len(chosen) == 2 and (chosen[0].label == "IR" or chosen[0].label != chosen[1].label):
                repeats = chosen
                break
    if len(repeats) != 2:
        return []
    first, second = sorted(repeats, key=lambda r: r.interval(length)[0])
    if first.touches(second.pieces, length):
        return []
    return _regions_between(first.interval(length), second.interval(length), first.label, second.label, length)


def _regions_between(a: tuple[int, int, int], b: tuple[int, int, int], a_label: str, b_label: str,
                     length: int) -> list[Region]:
    """The four regions made by two inverted repeats `a` and `b` ((start, end, length), `a` starting first): the
    repeats, named as given when they are IRa and IRb, else by convention (IRb follows the LSC), and the
    single-copy regions as the gaps between them, the larger one being the LSC. Nothing when a gap is empty."""
    a_start, a_end, a_len = a
    b_start, b_end, b_len = b
    gap_a = (a_end % length + 1, (b_start - 2) % length + 1, (b_start - a_end - 1) % length)  # Between them
    gap_b = (b_end % length + 1, (a_start - 2) % length + 1, (a_start - b_end - 1) % length)  # Around the origin
    len_a, len_b = gap_a[2], gap_b[2]
    if len_a <= 0 or len_b <= 0 or a_len + b_len + len_a + len_b != length:
        return []
    if {a_label, b_label} != {"IRa", "IRb"}:
        a_label, b_label = ("IRa", "IRb") if len_a >= len_b else ("IRb", "IRa")
    lsc_first = len_a >= len_b
    regions = [Region(a_label, a_start, a_end, a_len), Region("LSC" if lsc_first else "SSC", *gap_a),
               Region(b_label, b_start, b_end, b_len), Region("SSC" if lsc_first else "LSC", *gap_b)]
    return sorted(regions, key=lambda r: r.start)


# ---------------------------------------------------------------------------------------------------------------
# The inverted repeat of a plastome, detected in its sequence
# ---------------------------------------------------------------------------------------------------------------

IR_K = 32  # Seeds: k-mers of one strand matching k-mers of the other
IR_STEP = 4  # Every IR_STEP-th k-mer of the sequence is indexed; every k-mer of the other strand is looked up
IR_BAND = 500  # Seeds within this many antidiagonals of the chain's belong to the same repeat (indels shift them)
IR_MAX_GAP = 2000  # A longer stretch of a copy without a seed ends the repeat
IR_MAX_INDEL = 3000  # The chain goes on past an indel of up to this many bases when the seeds resume beyond it
IR_LOOKAHEAD = 12  # At the ends, a mismatch is passed when this many bases match after it
IR_FEW = 3  # Equal stretches with more mismatches than this are aligned (compensating indels look like mismatches)
IR_MARGIN = 16  # Stretches are aligned within a band of their length difference plus this many bases either side
IR_MAX_CELLS = 400_000  # A larger alignment (cells) is not done: the stretch counts as all different
MIN_DETECTED_REPEAT = 5000  # Each copy of a detected repeat must be this long (plastid IRs are 10-30 kb)
MIN_IDENTITY = 0.99  # ... and the copies this identical
IR_PENALTY = 20  # Trimming: a base of the chain scores +1, a difference -20 (a flank under 95% identical, < 0)
MAX_DETECTION_LENGTH = 2_000_000  # Longer sequences (bacteria) are not searched
AGREE_FRACTION = 0.8  # A detected copy lying this much inside an annotated copy confirms the annotation
PLASTOME_MIN_IR_FRACTION = 0.05  # A plastome layout: the two repeats make at least this much of the sequence ...
PLASTOME_MAX_SINGLE_COPY = 200_000  # ... and the larger single-copy region is at most this long
PLASTID_GENOMES = {"plastid", "chloroplast", "apicoplast", "cyanelle", "chromoplast", "leucoplast", "proplastid"}
_IR_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_NON_ACGT = re.compile("[^ACGT]+")


@dataclass
class InvertedRepeat:
    """Two copies of an inverted repeat found in a sequence: 1-based inclusive coordinates (the end smaller than
    the start when a copy spans the origin), the copies' lengths and the differences between them: mismatches,
    indels and runs of bases other than ACGT, each counted once whatever its length."""

    first: tuple[int, int]  # The copy starting first in the sequence
    second: tuple[int, int]
    lengths: tuple[int, int]
    differences: int

    @property
    def identity(self) -> float:
        return 1 - self.differences / max(self.lengths)

    def text(self) -> str:
        """`two copies of 25,593 bp, 100% identical` (`copies of 25,593 and 25,591 bp, 99.96% identical` when they
        differ)."""
        a, b = self.lengths
        size = f"two copies of {a:,} bp" if a == b else f"copies of {a:,} and {b:,} bp"
        pct = "100%" if not self.differences else f"{min(self.identity * 100, 99.99):.2f}%".replace(".00%", "%")
        return f"{size}, {pct} identical"


def _repeat_seeds(seq: str, k: int, step: int) -> list[tuple[int, int]]:
    """Pairs (i, j), i < j, of 0-based positions whose k-mers are reverse complements of each other: every
    step-th k-mer of the sequence is indexed (once: a k-mer at several indexed positions, or with a base other
    than ACGT, is left out), and every k-mer of the reverse complement is looked up."""
    n = len(seq)
    masked = {i for m in _NON_ACGT.finditer(seq) for i in range(max(0, m.start() - k + 1), m.end())}
    index: dict[str, int] = {}
    for i in range(0, n - k + 1, step):
        if i not in masked:
            kmer = seq[i:i + k]
            index[kmer] = -1 if kmer in index else i
    rc = seq.translate(_IR_COMPLEMENT)[::-1]
    get = index.get
    seeds = set()
    for p in range(n - k + 1):
        i = get(rc[p:p + k])
        if i is not None and i >= 0:
            j = n - p - k  # The k-mer at j, read on the other strand, is the k-mer at i
            if i != j:
                seeds.add((i, j) if i < j else (j, i))
    return sorted(seeds)


def _band(seeds: list[tuple[int, int]], centre: int) -> list[tuple[int, int]]:
    """The seeds within IR_BAND of the antidiagonal `centre` (i + j), in order."""
    return [(i, j) for i, j in seeds if abs(i + j - centre) <= IR_BAND]


def _walk(band: list[tuple[int, int]], k: int) -> list[list[tuple[int, int]]]:
    """Chain the seeds of one band in order (i increasing, j decreasing), a gap of more than IR_MAX_GAP without
    a seed breaking the chain. A seed on another antidiagonal is only taken when the current one has no seed
    within IR_MAX_GAP after it (a short duplication inside the repeat seeds a parallel antidiagonal; a real
    indel ends the current one)."""
    by_antidiagonal: dict[int, list[int]] = {}
    for i, j in band:
        by_antidiagonal.setdefault(i + j, []).append(i)
    chains: list[list[tuple[int, int]]] = []
    chain: list[tuple[int, int]] = []
    for i, j in band:
        if chain:
            pi, pj = chain[-1]
            if j >= pj or i <= pi:
                continue
            if i - pi > IR_MAX_GAP:
                chains.append(chain)
                chain = []
            elif i + j != pi + pj:
                same = by_antidiagonal[pi + pj]
                nxt = bisect_right(same, pi)
                if nxt < len(same) and same[nxt] - pi <= IR_MAX_GAP:
                    continue
        chain.append((i, j))
    chains.append(chain)
    return chains


def _links(g1: int, g2: int) -> bool:
    """Whether two seeds `g1` bases apart in one copy and `g2` in the other can follow each other in a chain: at
    most IR_MAX_GAP bases of the shorter stretch without a seed, and an indel of at most IR_MAX_INDEL."""
    return min(g1, g2) <= IR_MAX_GAP and abs(g1 - g2) <= IR_MAX_INDEL


def _continue(seeds: list[tuple[int, int]], chain: list[tuple[int, int]], k: int) -> list[tuple[int, int]]:
    """Extend a chain forwards past indels larger than IR_BAND: when a seed beyond the chain's last one links to
    it (_links) on another antidiagonal (the one with the smallest indel, then the nearest), the chain goes on
    along that seed's band, and so on."""
    while True:
        pi, pj = chain[-1]
        lo = bisect_right(seeds, (pi, pj))
        hi = bisect_left(seeds, (pi + k + IR_MAX_GAP + IR_MAX_INDEL + 1, -1))
        found = None
        for i, j in seeds[lo:hi]:
            if j >= pj:
                continue
            g1, g2 = i - pi - k, pj - j - k
            if _links(g1, g2) and (found is None or (abs(g1 - g2), i) < found[0]):
                found = ((abs(g1 - g2), i), (i, j))
        if found is None:
            return chain
        i, j = found[1]
        chain = chain + _walk([(i, j), *(s for s in _band(seeds, i + j) if s[0] > i)], k)[0]


def _repeat_chain(seeds: list[tuple[int, int]], k: int) -> list[tuple[int, int]]:
    """The chain of seeds of the largest inverted repeat: the longest chain (_walk) through the antidiagonal
    (i + j constant) with the most seeds and those within IR_BAND of it, continued at both ends past larger
    indels (_continue)."""
    if not seeds:
        return []
    counts = Counter(i + j for i, j in seeds)
    best = max(counts, key=lambda c: (counts[c], -c))
    chain = max(_walk(_band(seeds, best), k), key=lambda ch: ch[-1][0] + k - ch[0][0] if ch else 0)
    if not chain:
        return []
    chain = _continue(seeds, chain, k)
    mirrored = sorted((-i, -j) for i, j in seeds)  # Backwards from the chain's start is forwards here
    chain = _continue(mirrored, [(-i, -j) for i, j in reversed(chain)], k)
    return [(-i, -j) for i, j in reversed(chain)]


def _banded_edit_distance(a: str, b: str, band: int) -> int:
    """The edit distance of `a` to `b` (len(a) <= len(b); a base other than ACGT matches anything) through the
    cells within `band` of the diagonals running from the first bases to the last."""
    la, lb = len(a), len(b)
    delta = lb - la
    inf = la + lb + 1
    wild = bool(_NON_ACGT.search(a + b))
    prev = list(range(min(delta + band, lb) + 1))  # Row 0: y insertions
    plo = 0
    for x in range(1, la + 1):
        lo, hi = max(0, x - band), min(lb, x + delta + band)
        up = prev[lo - plo:]  # prev at (x - 1, y) for y from lo; beyond its end: out of the band
        if len(up) < hi - lo + 1:
            up.append(inf)
        diag = prev[lo - 1 - plo] if lo - 1 >= plo else inf  # prev at (x - 1, lo - 1)
        ca = a[x - 1]
        wild_a = wild and ca not in "ACGT"
        if lo == 0:
            cur, left, start = [x], x, 1
            diag = up[0]
        else:
            cur, left, start = [], inf, lo
        append = cur.append
        for y in range(start, hi + 1):
            u = up[y - lo]
            cb = b[y - 1]
            v = diag if (ca == cb or wild_a or (wild and cb not in "ACGT")) else diag + 1
            if u + 1 < v:
                v = u + 1
            if left + 1 < v:
                v = left + 1
            append(v)
            left, diag = v, u
        prev, plo = cur, lo
    return prev[lb - plo]


def _matches(x: str, y: str) -> bool:
    return x == y or x not in "ACGT" or y not in "ACGT"


def _ambiguous_runs(a: str, b: str) -> int:
    """The runs of columns of two aligned stretches (same length) differing by a base other than ACGT: an N run
    in one copy counts once, the same ambiguity code in both copies not at all."""
    runs, inside = 0, False
    for x, y in zip(a, b):
        odd = x != y and (x not in "ACGT" or y not in "ACGT")
        if odd and not inside:
            runs += 1
        inside = odd
    return runs


def _stretch_differences(a: str, b: str) -> int:
    """The differences between the bases of each copy between two seeds (`b` read on the first copy's strand):
    mismatches, indels and runs of bases other than ACGT (which match anything), each counted once. After their
    common ends are removed, stretches of equal length are compared base by base unless that gives more than
    IR_FEW mismatches (compensating indels look like a stretch of mismatches): then, like stretches of different
    lengths, they are aligned within a band of their length difference plus IR_MARGIN; an alignment of more than
    IR_MAX_CELLS cells is not done and the stretch counts as all different."""
    ambiguous = bool(_NON_ACGT.search(a + b))
    p = 0
    while p < len(a) and p < len(b) and _matches(a[p], b[p]):
        p += 1
    s = 0
    while s < len(a) - p and s < len(b) - p and _matches(a[-1 - s], b[-1 - s]):
        s += 1
    runs = _ambiguous_runs(a[:p], b[:p]) + _ambiguous_runs(a[len(a) - s:], b[len(b) - s:]) if ambiguous else 0
    a, b = a[p:len(a) - s], b[p:len(b) - s]
    if len(a) > len(b):
        a, b = b, a
    delta = len(b) - len(a)
    if not delta:
        if ambiguous:
            runs += _ambiguous_runs(a, b)
        mismatches = sum(not _matches(x, y) for x, y in zip(a, b))
        if mismatches <= IR_FEW:
            return mismatches + runs
    elif ambiguous:
        runs += len(_NON_ACGT.findall(a)) + len(_NON_ACGT.findall(b))
    if not a:
        return 1 + runs  # An indel
    if len(a) * (delta + 2 * IR_MARGIN + 1) > IR_MAX_CELLS:
        return len(a) + 1 + runs
    return _banded_edit_distance(a, b, IR_MARGIN) - max(0, delta - 1) + runs


def _link_differences(seq: str, first: tuple[int, int], second: tuple[int, int], k: int) -> int:
    """The differences between the copies between two consecutive seeds of a chain: those of the stretches
    between them (_stretch_differences), or one indel when the seeds overlap by different amounts."""
    (i1, j1), (i2, j2) = first, second
    g1, g2 = i2 - i1 - k, j1 - j2 - k  # The bases of each copy between the two seeds
    if g1 <= 0 or g2 <= 0:
        return int(g1 != g2)
    return _stretch_differences(seq[i1 + k:i2], reverse_complement(seq[j2 + k:j1]))


def _best_subchain(chain: list[tuple[int, int]], differences: list[int], k: int) -> tuple[int, int, int]:
    """The part of a chain scoring most, a base of the first copy counting +1 and a difference -IR_PENALTY (a
    diverged flank scores less than nothing): the indexes of its first and last seeds and its differences."""
    best = (k, 0, 0, 0)  # Score, first seed, last seed, differences
    running, start, count = k, 0, 0
    for idx, ((i1, _), (i2, _)) in enumerate(zip(chain, chain[1:])):
        gain = i2 - i1 - IR_PENALTY * differences[idx]
        if running + gain < k:  # Starting afresh at the next seed scores more
            running, start, count = k, idx + 1, 0
        else:
            running, count = running + gain, count + differences[idx]
        if running >= best[0]:  # The longer chain on a tie
            best = (running, start, idx + 1, count)
    return best[1], best[2], best[3]


def _shift_interval(interval: tuple[int, int], offset: int, n: int) -> tuple[int, int]:
    s, e = interval
    return (s - 1 + offset) % n + 1, (e - 1 + offset) % n + 1


def _may_continue(seq: str, a: int, b_end: int, k: int) -> bool:
    """Whether copies reaching the origin (the first starting at the first base, or the second ending at the
    last) may go on beyond it after an indel: a k-mer of the bases after the second copy matches, read on the
    other strand, one of the bases before the first (IR_MAX_INDEL + k bases each, fewer when the copies would
    meet)."""
    n = len(seq)
    w = min(IR_MAX_INDEL + k, (a + n - b_end) // 2)
    if w < k:
        return False
    after = (seq[b_end:] + seq)[:w]
    before = reverse_complement((seq + seq[:a])[a - w + n:a + n])
    kmers = {after[p:p + k] for p in range(w - k + 1)}
    return any(before[p:p + k] in kmers for p in range(w - k + 1))


def _same_pair(found: InvertedRepeat, repeat: InvertedRepeat, n: int) -> bool:
    """Whether two detections are of the same repeat: their copies overlap, pair by pair (in either pairing: a
    copy across the origin starts late in the sequence but covers its first base)."""
    f = [_pieces(*c, n) for c in (found.first, found.second)]
    r = [_pieces(*c, n) for c in (repeat.first, repeat.second)]
    return any(_overlap(f[0], r[x]) > 0 and _overlap(f[1], r[1 - x]) > 0 for x in (0, 1))


def detect_inverted_repeat(seq: str, _rotated: bool = False) -> InvertedRepeat | None:  # noqa: C901
    """The large inverted repeat of a sequence (circular: a copy may span the origin), or None: two copies of at
    least MIN_DETECTED_REPEAT bp each, at least MIN_IDENTITY identical, not overlapping (they may abut: the two
    halves of a palindrome are two copies). Seeds (k-mers matching k-mers of the other strand, see _repeat_seeds)
    are chained (_repeat_chain), the differences between the copies are counted between the seeds (mismatches,
    indels and runs of N, each once: _stretch_differences), the chain is trimmed to its part scoring most
    (_best_subchain, which drops diverged flanks) and extended base by base at both ends (through a mismatch
    followed by IR_LOOKAHEAD matching bases). A repeat found across the origin of the sequence, or ending at it
    with matching k-mers beyond (_may_continue), is searched again in the sequence rotated to start between the
    copies, where the whole of each copy is seeded; that result is taken when it is the same pair of copies.
    Sequences shorter than two copies or longer than MAX_DETECTION_LENGTH are not searched. Linear in the length:
    about a tenth of a second for a plastome, a second for 2 Mb."""
    seq = seq.upper()
    n = len(seq)
    k = IR_K
    if n < 2 * MIN_DETECTED_REPEAT or n > MAX_DETECTION_LENGTH:
        return None
    chain = _repeat_chain(_repeat_seeds(seq, k, IR_STEP), k)
    if not chain:
        return None
    per_link = [_link_differences(seq, s, t, k) for s, t in zip(chain, chain[1:])]
    start, end, differences = _best_subchain(chain, per_link, k)
    chain = chain[start:end + 1]
    a, a_end = chain[0][0], chain[-1][0] + k  # The first copy, [a, a_end)
    b, b_end = chain[-1][1], chain[0][1] + k  # The second, [b, b_end)
    if a_end > b:  # The innermost seeds overlap (by less than k): the copies abut
        a_end = b = (a_end + b) // 2
    look = IR_LOOKAHEAD

    def matches(x: int, y: int) -> bool:
        p, q = seq[x % n], seq[y % n]
        return p in "ACGT" and p == q.translate(_IR_COMPLEMENT)

    # Outwards (the first copy leftwards, the second rightwards, around the origin when needed) and inwards; the
    # copies must not meet
    while a + n - b_end >= 2:
        if matches(a - 1, b_end):
            a, b_end = a - 1, b_end + 1
        elif a + n - b_end > 2 * (look + 1) and all(matches(a - 1 - t, b_end + t) for t in range(1, look + 1)):
            a, b_end, differences = a - 1, b_end + 1, differences + 1
        else:
            break
    while b - a_end >= 2:
        if matches(a_end, b - 1):
            a_end, b = a_end + 1, b - 1
        elif b - a_end > 2 * (look + 1) and all(matches(a_end + t, b - 1 - t) for t in range(1, look + 1)):
            a_end, b, differences = a_end + 1, b - 1, differences + 1
        else:
            break
    repeat = InvertedRepeat((a % n + 1, (a_end - 1) % n + 1), (b % n + 1, (b_end - 1) % n + 1),
                            (a_end - a, b_end - b), differences)
    if not _rotated and (a < 0 or b_end > n or ((a == 0 or b_end == n) and _may_continue(seq, a, b_end, k))):
        offset = (a_end + b) // 2 % n
        found = detect_inverted_repeat(seq[offset:] + seq[:offset], _rotated=True)
        if found is not None:
            copies = sorted(zip((_shift_interval(c, offset, n) for c in (found.first, found.second)),
                                found.lengths))
            found = InvertedRepeat(copies[0][0], copies[1][0], (copies[0][1], copies[1][1]), found.differences)
            if _same_pair(found, repeat, n):
                return found
    if min(repeat.lengths) < MIN_DETECTED_REPEAT or repeat.identity < MIN_IDENTITY:
        return None
    return repeat


def regions_from_repeat(repeat: InvertedRepeat, length: int) -> list[Region]:
    """The LSC/IRb/SSC/IRa regions of a sequence from its detected inverted repeat, named by convention (IRb
    follows the LSC, the larger single-copy region); none when the copies abut."""
    copies = sorted((repeat.first, repeat.second), key=lambda c: c[0])
    intervals = [(s, e, repeat.lengths[(repeat.first, repeat.second).index((s, e))]) for s, e in copies]
    return _regions_between(intervals[0], intervals[1], "IR", "IR", length)


@dataclass
class RegionBand:
    """The LSC/IRb/SSC/IRa regions of one sequence and where they come from: the annotated inverted repeats
    (`annotation`) or the inverted repeat detected in the sequence (`sequence`); or no regions (`none`) with
    the repeat found and a note saying why it gives none."""

    regions: list[Region]
    source: str
    repeat: InvertedRepeat | None = None
    note: str = ""

    def text(self) -> str:
        """`the annotated inverted repeats` or `the inverted repeat detected in the reference sequence (two copies
        of 25,593 bp, 100% identical)`."""
        if self.source == "annotation":
            return "the annotated inverted repeats"
        detail = f" ({self.repeat.text()})" if self.repeat else ""
        return f"the inverted repeat detected in the reference sequence{detail}"

    def record(self) -> dict[str, object]:
        """For run_info.json."""
        item: dict[str, object] = {
            "source": self.source,
            "regions": [{"name": r.name, "start": r.start, "end": r.end, "length": r.length}
                        for r in self.regions]}
        if self.repeat:
            item["repeat"] = {"copies": [list(self.repeat.first), list(self.repeat.second)],
                              "lengths": list(self.repeat.lengths), "differences": self.repeat.differences,
                              "identity": round(self.repeat.identity, 6)}
        if self.note:
            item["note"] = self.note
        return item


def _pieces(start: int, end: int, length: int) -> list[tuple[int, int]]:
    return [(start, end)] if start <= end else [(start, length), (1, end)]


def _overlap(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    return sum(max(0, min(e1, e2) - max(s1, s2) + 1) for s1, e1 in a for s2, e2 in b)


def _agree(annotated: list[Region], repeat: InvertedRepeat, length: int) -> bool:
    """Whether the detected repeat confirms the annotated one: each detected copy lies at least AGREE_FRACTION
    inside an annotated copy (a detection cut short by a large insertion or a run of N in one copy confirms the
    annotation; copies found elsewhere, as when the annotation names the single-copy regions as the repeats,
    contradict it)."""
    pieces = [_pieces(r.start, r.end, length) for r in annotated if r.name.startswith("IR")]
    for copy, size in zip((repeat.first, repeat.second), repeat.lengths):
        detected = _pieces(*copy, length)
        if not any(_overlap(detected, p) >= AGREE_FRACTION * size for p in pieces):
            return False
    return True


def _not_a_plastid(features: list[RawFeature]) -> str:
    """Why the annotation says the sequence is not a plastid genome (`the source feature says
    organelle=mitochondrion`, `the region feature says genome=chromosome`), or '' when it says it is one, or
    nothing: a GenBank source feature's /organelle, or an NCBI GFF3 region's genome= attribute."""
    for f in features:
        if f.type == "source" and f.qualifiers.get("organelle"):
            value = f.qualifiers["organelle"].split("; ")[0]
            kind = value.lower().partition(":")[2] or value.lower()
            return "" if kind in PLASTID_GENOMES else f"the source feature says organelle={value}"
        if f.type == "region" and f.qualifiers.get("genome"):
            value = f.qualifiers["genome"].split("; ")[0]
            return "" if value.lower() in PLASTID_GENOMES else f"the region feature says genome={value}"
    return ""


def _layout_problem(regions: list[Region], length: int, not_plastid: str = "") -> str:
    """Why four regions are not the quadripartite layout of a plastome ('' when they are): the sequence is
    said not to be a plastid genome, the repeats make less than PLASTOME_MIN_IR_FRACTION of the sequence (the
    inverted rRNA operons of a bacterium, a small repeat of a plant mitochondrion), or the larger single-copy
    region is longer than PLASTOME_MAX_SINGLE_COPY."""
    if not_plastid:
        return not_plastid
    repeats = sum(r.length for r in regions if r.name.startswith("IR"))
    single = max(r.length for r in regions if not r.name.startswith("IR"))
    if repeats < PLASTOME_MIN_IR_FRACTION * length:
        return (f"the repeats are {repeats / length:.1%} of the sequence, less than "
                f"{PLASTOME_MIN_IR_FRACTION:.0%}")
    if single > PLASTOME_MAX_SINGLE_COPY:
        return (f"the larger single-copy region is {single / 1000:,.0f} kb, more than "
                f"{PLASTOME_MAX_SINGLE_COPY // 1000} kb")
    return ""


def find_regions(features: list[RawFeature], length: int, seq: str | None = None, name: str = "",  # noqa: C901
                 warnings: list[str] | None = None) -> RegionBand | None:
    """The regions of a sequence: from its annotated inverted repeats (derive_regions) when they make a
    plastome layout (_layout_problem) that the inverted repeat found in the sequence confirms (_agree), or when
    the sequence is not given or has no detectable repeat; else from the detected repeat (regions_from_repeat),
    with a warning (in `warnings`) when the annotation named inverted repeats. A repeat, annotated or detected,
    that is not a plastome layout gives no regions: a RegionBand without regions (`none`) keeps the detected
    repeat and a note saying why, for run_info.json. None when nothing was found."""
    annotated = derive_regions(features, length) if features else []
    named = bool(features) and bool(_inverted_repeats(features, length))
    not_plastid = _not_a_plastid(features) if features else ""
    annotated_problem = _layout_problem(annotated, length, not_plastid) if annotated else ""
    repeat = detect_inverted_repeat(seq) if seq is not None else None
    detected = regions_from_repeat(repeat, length) if repeat else []
    detected_problem = "" if not repeat else ("the copies abut" if not detected
                                              else _layout_problem(detected, length, not_plastid))
    copies = ", ".join(f"{r.start:,}–{r.end:,}" for r in annotated if r.name.startswith("IR"))
    if annotated and not annotated_problem and (repeat is None or detected_problem
                                                or _agree(annotated, repeat, length)):
        return RegionBand(annotated, "annotation")
    found = "" if repeat is None else (f"{repeat.first[0]:,}–{repeat.first[1]:,} and {repeat.second[0]:,}–"
                                       f"{repeat.second[1]:,}")
    if repeat is None or detected_problem:
        if warnings is not None and annotated:  # Then annotated_problem: a plausible annotation was returned
            warnings.append(f"{name}: the annotated inverted repeats ({copies}) are not a plastome layout "
                            f"({annotated_problem}); no regions from them")
        if repeat is None:
            return RegionBand([], "none", None, f"annotated inverted repeats but not a plastome layout "
                                                f"({annotated_problem})") if annotated else None
        return RegionBand([], "none", repeat, f"inverted repeat found but not a plastome layout "
                                              f"({detected_problem})")
    if warnings is not None and annotated and annotated_problem:
        warnings.append(f"{name}: the annotated inverted repeats ({copies}) are not a plastome layout "
                        f"({annotated_problem}); the regions follow the inverted repeat found in the sequence "
                        f"({found})")
    elif warnings is not None and annotated:
        warnings.append(f"{name}: the annotated inverted repeats ({copies}) are not the inverted repeat found in "
                        f"the sequence ({found}); the regions follow the sequence")
    elif warnings is not None and named:
        warnings.append(f"{name}: the annotated inverted repeats do not give a plastome layout; the regions "
                        f"follow the inverted repeat found in the sequence ({found})")
    return RegionBand(detected, "sequence", repeat)


# ---------------------------------------------------------------------------------------------------------------
# The annotation of a reference
# ---------------------------------------------------------------------------------------------------------------

class SequenceAnnotation:
    """The genes and regions of one reference sequence, indexed by position."""

    def __init__(self, name: str, length: int, genes: list[Gene], regions: list[Region] | RegionBand | None,
                 circular: bool):
        self.name, self.length, self.circular = name, length, circular
        self.genes = sorted(genes, key=lambda g: (g.start, g.end))
        if isinstance(regions, RegionBand):
            self.band: RegionBand | None = regions
        else:
            self.band = RegionBand(regions, "annotation") if regions else None
        self.regions = self.band.regions if self.band else []
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


def load_annotation(path: Path, sequences: list[tuple[str, int]], default_table: int = DEFAULT_TABLE,
                    seqs: dict[str, str] | None = None) -> Annotation:
    """Read a GenBank or GFF3 annotation and match it to the reference sequences (name, length). An unreadable
    file or an unknown format is a BaconError; annotated sequences that match no reference sequence, and
    features beyond the end of their sequence, are left out with a warning (in `warnings`). A single annotated
    sequence is taken for a single reference sequence of the same length, whatever its name. With the reference
    sequences themselves (`seqs`, by name), the regions of an annotated sequence without annotated inverted
    repeats come from the inverted repeat detected in its sequence (find_regions)."""
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
        features, lengths, circular = read_gff3(path, dict(sequences))
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
    for old, new in rename.items():
        if lengths.get(old) and lengths[old] != reference[new]:
            warnings.append(f"{path.name}: the annotation of {new} is {lengths[old]:,} bp, the reference "
                            f"{reference[new]:,} bp: is it the annotation of this reference?")
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
            band = find_regions(kept[name], length, (seqs or {}).get(name), name, warnings)
            annotation.sequences[name] = SequenceAnnotation(
                name, length, builder(kept[name], default_table, length), band, is_circular)
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
# Differences from the standard code and the start codons of the known tables, as NCBI defines them
# (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi): 1 standard (with ATG only), 2 vertebrate, 3 yeast,
# 4 mold/protozoan, 5 invertebrate, 9 echinoderm and flatworm, 13 ascidian and 14 alternative flatworm
# mitochondrial, 11 bacterial, archaeal and plant plastid.
TABLES: dict[int, tuple[dict[str, str], set[str]]] = {
    1: ({}, {"ATG"}),
    2: ({"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W"}, {"ATT", "ATC", "ATA", "ATG", "GTG"}),
    3: ({"ATA": "M", "CTT": "T", "CTC": "T", "CTA": "T", "CTG": "T", "TGA": "W"}, {"ATA", "ATG", "GTG"}),
    4: ({"TGA": "W"}, {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA", "TTA"}),
    5: ({"AGA": "S", "AGG": "S", "ATA": "M", "TGA": "W"}, {"TTG", "ATT", "ATC", "ATA", "ATG", "GTG"}),
    9: ({"AAA": "N", "AGA": "S", "AGG": "S", "TGA": "W"}, {"ATG", "GTG"}),
    11: ({}, {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}),
    13: ({"AGA": "G", "AGG": "G", "ATA": "M", "TGA": "W"}, {"TTG", "ATA", "ATG", "GTG"}),
    14: ({"AAA": "N", "AGA": "S", "AGG": "S", "TAA": "Y", "TGA": "W"}, {"ATG"}),
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
    (a partial CDS), contains N or has a translational exception (/transl_except: selenocysteine, an edited
    codon), or the VCF's reference base does not match the reference."""
    if len(ref) != 1 or len(alt) != 1 or alt.upper() not in "ACGT":
        return None
    found = cds.index(pos)
    if found is None or found[0] < 0 or cds.excepted(pos):
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
    """The region, genes, context and effects of a SNP. The same effect through two coding sequences of a gene
    (an exon shared by both products of a trans-spliced rps12) is given once."""
    genes = seq_ann.genes_at(pos)
    region = seq_ann.region_at(pos)
    effects: list[Effect] = []
    if genes:
        contexts = list(dict.fromkeys(g.context(pos) for g in genes))
        seen: set[tuple[str, str, str, str, str]] = set()
        for g in genes:
            if g.kind == "pseudogene":
                continue
            for cds in g.cds:
                for alt in alts:
                    effect = cds_effect(cds, g.name, pos, ref, alt, seq)
                    if effect is None:
                        continue
                    key = (effect.gene, effect.alt, effect.codons, effect.change, effect.kind)
                    if key not in seen:
                        seen.add(key)
                        effects.append(effect)
        return SnpAnnotation(region, genes, " / ".join(contexts), effects)
    before, after = seq_ann.neighbours(pos)
    if before is not None and after is not None:
        context = f"intergenic between {before.label} and {after.label}"
    elif after is not None:
        context = f"intergenic before {after.label}"
    elif before is not None:
        context = f"intergenic after {before.label}"
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
