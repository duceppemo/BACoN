"""Finding the samples to process: an input file, an input folder, or a sample sheet."""

from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass
from pathlib import Path

from bacon import BaconError
from bacon.seqio import sniff_format, split_extension

log = logging.getLogger(__name__)

# Folders that MinKNOW/Dorado write next to the barcode folders and that are not samples.
SKIPPED_FOLDERS = {"unclassified", "mixed"}
VALID_NAME = re.compile(r"^[A-Za-z0-9._+-]+$")


@dataclass
class Sample:
    name: str
    files: list[Path]
    fmt: str = "fastq"  # 'fastq' or 'fasta'


def _sequence_files(folder: Path) -> list[Path]:
    return sorted(p for p in folder.rglob("*")
                  if p.is_file() and not p.name.startswith(".") and split_extension(p.name))


def _check_name(name: str, origin: object) -> None:
    if not VALID_NAME.match(name):
        raise BaconError(f"Invalid sample name {name!r} (from {origin}): use only letters, digits and . _ + -")
    if name.lower() in {"reference", "all_assemblies", "assembly_graphs"}:
        raise BaconError(f"Sample name {name!r} (from {origin}) is reserved; rename the file or folder")


def _make_sample(name: str, files: list[Path], origin: object) -> Sample | None:
    _check_name(name, origin)
    formats = {f: sniff_format(f) for f in files}
    empty = [f for f, fmt in formats.items() if fmt is None]
    for f in empty:
        log.warning("Skipping empty file %s", f)
    files = [f for f in files if formats[f] is not None]
    if not files:
        log.warning("Skipping sample %s: no reads", name)
        return None
    kinds = {formats[f] for f in files}
    if len(kinds) > 1:
        raise BaconError(f"Sample {name} mixes fasta and fastq files: {', '.join(map(str, files))}")
    return Sample(name, files, kinds.pop())


def discover(path: Path) -> list[Sample]:
    """Samples from a single file, or from a folder.

    In a folder, each sequence file is one sample named after the file; each subfolder is one sample named after
    the subfolder, with all the sequence files it contains (for example MinKNOW's fastq_pass/barcode01/).
    """
    if not path.exists():
        raise BaconError(f"Input not found: {path}")
    candidates: list[tuple[str, list[Path], Path]] = []
    if path.is_file():
        parts = split_extension(path.name)
        if not parts:
            raise BaconError(f"Not a fasta/fastq file (by extension): {path}")
        candidates.append((parts[0], [path], path))
    else:
        for entry in sorted(path.iterdir()):
            if entry.name.startswith("."):
                continue
            if entry.is_dir():
                files = _sequence_files(entry)
                if not files:
                    continue
                if entry.name in SKIPPED_FOLDERS:
                    log.info("Skipping folder %s (reads without a barcode)", entry)
                    continue
                candidates.append((entry.name, files, entry))
            elif entry.is_file() and (parts := split_extension(entry.name)):
                candidates.append((parts[0], [entry], entry))
    return _finalize(candidates, path)


def read_sample_sheet(path: Path) -> list[Sample]:
    """Samples from a TSV or CSV file with the columns 'sample' and 'file' (files separated by ';').

    Several rows may share a sample name; their files are merged. Relative paths are resolved from the
    sheet's folder.
    """
    if not path.is_file():
        raise BaconError(f"Sample sheet not found: {path}")
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    lines = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise BaconError(f"Sample sheet is empty: {path}")
    delimiter = "\t" if "\t" in lines[0] else ","
    reader = csv.DictReader(lines, delimiter=delimiter)
    columns = {c.strip().lower(): c for c in reader.fieldnames or []}
    if "sample" not in columns or "file" not in columns:
        raise BaconError(f"Sample sheet {path} needs the columns 'sample' and 'file' (found: {reader.fieldnames})")
    grouped: dict[str, list[Path]] = {}
    for row in reader:
        name = (row.get(columns["sample"]) or "").strip()
        files = [f.strip() for f in (row.get(columns["file"]) or "").split(";") if f.strip()]
        if not name or not files:
            raise BaconError(f"Sample sheet {path}: row without a sample name or a file: {row}")
        for f in files:
            file = Path(f).expanduser()
            file = file if file.is_absolute() else path.parent / file
            if not file.is_file():
                raise BaconError(f"Sample sheet {path}: file not found for {name}: {file}")
            grouped.setdefault(name, []).append(file)
    candidates = [(name, files, path) for name, files in grouped.items()]
    return _finalize(candidates, path)


def _finalize(candidates: list[tuple[str, list[Path], object]], source: Path) -> list[Sample]:
    seen: dict[str, object] = {}
    for name, _, origin in candidates:
        if name in seen:
            raise BaconError(f"Two inputs give the same sample name {name!r}: {seen[name]} and {origin}")
        seen[name] = origin
    samples = [s for name, files, origin in candidates if (s := _make_sample(name, files, origin))]
    if not samples:
        raise BaconError(f"No fasta/fastq files with reads found in {source}")
    return samples
