"""Reading and writing fasta/fastq files, gzipped or not."""

from __future__ import annotations

import gzip
import io
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

from bacon import BaconError

try:  # Optional: python-isal decompresses gzip several times faster
    from isal import igzip as _gzip
except ImportError:  # pragma: no cover - depends on the environment
    _gzip = gzip

FASTQ_EXTENSIONS = (".fastq", ".fq")
FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".fas")


def split_extension(filename: str) -> tuple[str, str] | None:
    """Return (stem, 'fastq' or 'fasta') for a sequence file name, or None if it is not one."""
    lower = filename.lower()
    for fmt, extensions in (("fastq", FASTQ_EXTENSIONS), ("fasta", FASTA_EXTENSIONS)):
        for ext in extensions:
            for suffix in (ext + ".gz", ext):
                if lower.endswith(suffix) and len(filename) > len(suffix):
                    return filename[: -len(suffix)], fmt
    return None


def open_text(path: Path) -> io.TextIOBase:
    """Open a text file for reading, decompressing it if it is gzipped (detected from its content)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return _gzip.open(path, "rt", encoding="ascii", errors="replace")
    return open(path, encoding="ascii", errors="replace")


def open_write(path: Path) -> io.TextIOBase:
    """Open a text file for writing, gzipped if its name ends with .gz."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt", compresslevel=4, encoding="ascii")
    return open(path, "w", encoding="ascii")


@dataclass
class Record:
    header: str  # Header line without the leading '>' or '@'
    seq: str
    qual: str | None = None

    @property
    def name(self) -> str:
        return self.header.split(maxsplit=1)[0] if self.header else ""

    def format(self) -> str:
        if self.qual is None:
            return f">{self.header}\n{self.seq}\n"
        return f"@{self.header}\n{self.seq}\n+\n{self.qual}\n"


def sniff_format(path: Path) -> str | None:
    """'fastq', 'fasta', or None for an empty file, from the first non-blank character."""
    with open_text(path) as fh:
        for line in fh:
            if line.strip():
                if line[0] == "@":
                    return "fastq"
                if line[0] == ">":
                    return "fasta"
                raise BaconError(f"{path} is neither fasta nor fastq (first line: {line.strip()[:40]!r})")
    return None


def read_records(path: Path) -> Iterator[Record]:
    """Iterate over the records of a fasta (multi-line allowed) or fastq (4 lines per record) file."""
    fmt = sniff_format(path)
    if fmt is None:
        return
    with open_text(path) as fh:
        if fmt == "fasta":
            header, chunks = None, []
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if header is not None:
                        yield Record(header, "".join(chunks))
                    header, chunks = line[1:], []
                elif line and header is not None:
                    chunks.append(line)
            if header is not None:
                yield Record(header, "".join(chunks))
            return
        while True:
            header = fh.readline()
            if not header:
                return
            if not header.strip():
                continue
            seq, plus, qual = fh.readline().rstrip(), fh.readline(), fh.readline().rstrip()
            if not header.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
                raise BaconError(f"{path}: malformed or truncated fastq record {header.strip()[:60]!r}")
            yield Record(header[1:].rstrip(), seq, qual)


@dataclass
class ReadStats:
    """Number of reads, bases and (optionally) read lengths, for N50."""

    reads: int = 0
    bases: int = 0
    lengths: list[int] = field(default_factory=list, repr=False)

    def add(self, length: int, keep_length: bool = True) -> None:
        self.reads += 1
        self.bases += length
        if keep_length:
            self.lengths.append(length)

    @property
    def n50(self) -> int:
        return n50(self.lengths)


def n50(lengths: Iterable[int]) -> int:
    ordered = sorted(lengths, reverse=True)
    half, running = sum(ordered) / 2, 0
    for length in ordered:
        running += length
        if running >= half:
            return length
    return 0


def read_stats(paths: Iterable[Path]) -> ReadStats:
    stats = ReadStats()
    for path in paths:
        for rec in read_records(path):
            stats.add(len(rec.seq))
    return stats


def sequence_lengths(path: Path) -> list[tuple[str, int]]:
    return [(rec.name, len(rec.seq)) for rec in read_records(path)]


def check_reference(path: Path) -> list[tuple[str, int]]:
    """Validate a fasta reference and return its (name, length) pairs."""
    if not path.is_file():
        raise BaconError(f"Reference file not found: {path}")
    if sniff_format(path) != "fasta":
        raise BaconError(f"The reference must be a fasta file: {path}")
    lengths = sequence_lengths(path)
    if not lengths or sum(n for _, n in lengths) == 0:
        raise BaconError(f"The reference contains no sequence: {path}")
    names = [name for name, _ in lengths]
    if len(set(names)) != len(names):
        raise BaconError(f"The reference has duplicate sequence names: {path}")
    return lengths


def write_fasta(path: Path, records: Iterable[Record], width: int = 80) -> int:
    """Write records as plain fasta wrapped at `width`. Returns the number of records written."""
    count = 0
    with open(path, "w", encoding="ascii") as out:
        for rec in records:
            out.write(f">{rec.header}\n")
            for i in range(0, len(rec.seq), width):
                out.write(rec.seq[i: i + width] + "\n")
            count += 1
    return count


def extract_reads(inputs: Iterable[Path], names: set[str], output: Path) -> tuple[ReadStats, ReadStats]:
    """Copy the reads whose name is in `names` from `inputs` to `output`.

    Returns statistics for all input reads (without lengths) and for the extracted reads.
    """
    total, kept = ReadStats(), ReadStats()
    with open_write(output) as out:
        for path in inputs:
            for rec in read_records(path):
                total.add(len(rec.seq), keep_length=False)
                if rec.name in names:
                    kept.add(len(rec.seq))
                    out.write(rec.format())
    return total, kept


def concatenate(inputs: Iterable[Path], output: Path) -> None:
    with open_write(output) as out:
        for path in inputs:
            for rec in read_records(path):
                out.write(rec.format())
