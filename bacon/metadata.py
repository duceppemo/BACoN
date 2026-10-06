"""Sample metadata (`--metadata`, or the extra columns of a sample sheet): read, merged, and the column that
colours the report's figures. Metadata is for the report only: it is not part of any checkpoint."""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path

from bacon import BaconError

log = logging.getLogger(__name__)

MISSING_VALUES = {"", "NA", "na", "-"}  # Cells read as "no value"
MAX_COLOUR_VALUES = 8  # Distinct values a column may have to colour the figures (the report's categorical slots)
MAX_VALUE_LENGTH = 30  # Above this average length of the distinct values, a column is free text
FREE_TEXT_FROM = 10  # With this many samples having a value, more distinct values than half of them is free text
KEY_COLUMNS = ("sample", "file")  # Columns of a sample sheet that are not metadata
COPY_NAME = "metadata.tsv"  # The normalised copy in the output folder


@dataclass
class Metadata:
    columns: list[str]
    rows: dict[str, dict[str, str]]  # Sample -> column -> value ("" when missing)
    warnings: list[str] = field(default_factory=list)

    def value(self, sample: str, column: str) -> str:
        return self.rows.get(sample, {}).get(column, "")

    def values(self, column: str, samples: list[str] | None = None) -> list[str]:
        """The values of a column (missing ones as ""), for every sample or for `samples`."""
        names = list(self.rows) if samples is None else samples
        return [self.value(n, column) for n in names]


def read_table(path: Path, what: str) -> tuple[list[str], list[dict[str, str]]]:
    """Header and rows of a TSV or CSV file: UTF-8 with or without BOM, LF or CRLF, lines starting with '#' and
    blank lines ignored; a tab in the header makes it a TSV, otherwise a CSV. Names and values are stripped;
    columns with an empty name (a trailing separator) are dropped."""
    if not path.is_file():
        raise BaconError(f"{what} not found: {path}")
    try:
        text = path.read_text(encoding="utf-8-sig", errors="replace")
    except OSError as exc:
        raise BaconError(f"{what} cannot be read: {path} ({exc.strerror})") from None
    lines = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise BaconError(f"{what} is empty: {path}")
    delimiter = "\t" if "\t" in lines[0] else ","
    reader = csv.DictReader(lines, delimiter=delimiter)
    header = [c.strip() for c in reader.fieldnames or []]
    rows = []
    for raw in reader:
        row = {}
        for key, value in raw.items():
            if key is None:  # More cells than header names
                continue
            if isinstance(value, list):  # pragma: no cover - csv gives a list for the restkey only
                value = value[0] if value else ""
            row[key.strip()] = " ".join((value or "").split())  # Tabs or line breaks would break the copy
        rows.append(row)
    return [c for c in header if c], rows


def _clean(value: str) -> str:
    return "" if value in MISSING_VALUES else value


def _find(header: list[str], name: str) -> str | None:
    """The header name matching `name` without regard to case."""
    return next((c for c in header if c.lower() == name), None)


def read_metadata(path: Path) -> Metadata:
    """The metadata file: a 'sample' column (any case) and any other columns. A sample listed twice keeps its
    first row (with a warning); rows without a sample name are ignored."""
    header, rows = read_table(path, "Metadata file")
    key = _find(header, "sample")
    if key is None:
        raise BaconError(f"Metadata file {path} needs a 'sample' column (found: {', '.join(header) or 'nothing'})")
    columns = [c for c in header if c != key]
    if not columns:
        raise BaconError(f"Metadata file {path} has no column besides {key!r}")
    metadata = Metadata(columns, {})
    duplicates: list[str] = []
    for row in rows:
        name = row.get(key, "")
        if not name:
            continue
        if name in metadata.rows:
            duplicates.append(name)
            continue
        metadata.rows[name] = {c: _clean(row.get(c, "")) for c in columns}
    if duplicates:
        shown = ", ".join(dict.fromkeys(duplicates))
        metadata.warnings.append(f"Metadata file {path}: {len(set(duplicates))} sample(s) listed more than once, "
                                 f"the first row is used: {shown}")
    return metadata


def sheet_metadata(path: Path) -> Metadata | None:
    """The metadata in a sample sheet: its columns other than 'sample' and 'file'. A sample on several rows keeps
    the first value of each column; rows giving different values are reported."""
    header, rows = read_table(path, "Sample sheet")
    key = _find(header, "sample")
    columns = [c for c in header if c.lower() not in KEY_COLUMNS]
    if key is None or not columns:
        return None
    metadata = Metadata(columns, {})
    conflicts: list[str] = []
    for row in rows:
        name = row.get(key, "")
        if not name:
            continue
        values = metadata.rows.setdefault(name, {c: "" for c in columns})
        for c in columns:
            value = _clean(row.get(c, ""))
            if not value:
                continue
            if not values[c]:
                values[c] = value
            elif values[c] != value:
                conflicts.append(f"{name} {c} ({values[c]!r} kept, not {value!r})")
    if conflicts:
        metadata.warnings.append(f"Sample sheet {path}: rows of the same sample give different values; the first "
                                 f"one is used: {'; '.join(conflicts[:5])}{' …' if len(conflicts) > 5 else ''}")
    return metadata


def merge(first: Metadata | None, second: Metadata | None) -> Metadata | None:
    """One table from two: the columns of `first` (all its values), then the columns of `second` that `first`
    lacks."""
    if first is None or second is None:
        return first or second
    columns = first.columns + [c for c in second.columns if c not in first.columns]
    rows = {}
    for name in dict.fromkeys([*first.rows, *second.rows]):
        rows[name] = {c: (first.value(name, c) if c in first.columns else second.value(name, c)) for c in columns}
    return Metadata(columns, rows, first.warnings + second.warnings)


def restrict(metadata: Metadata, samples: list[str]) -> tuple[Metadata, list[str]]:
    """The metadata of the run's samples, in their order (samples without a row get blank cells), and the names
    of the rows that match no sample."""
    rows = {name: {c: metadata.value(name, c) for c in metadata.columns} for name in samples}
    unmatched = [name for name in metadata.rows if name not in rows]
    return Metadata(metadata.columns, rows, metadata.warnings), unmatched


def unusable_reason(values: list[str]) -> str | None:
    """Why a column cannot colour the figures, or None when it can: it needs a value, at most MAX_COLOUR_VALUES
    distinct ones, and must not look like free text (long values, or nearly one value per sample)."""
    present = [v for v in values if v]
    distinct = sorted(set(present))
    if not distinct:
        return "no value"
    if len(distinct) > MAX_COLOUR_VALUES:
        return f"{len(distinct)} distinct values (at most {MAX_COLOUR_VALUES} can be coloured)"
    if len(present) >= FREE_TEXT_FROM and len(distinct) > len(present) / 2:
        return f"{len(distinct)} distinct values among {len(present)} samples (free text)"
    if sum(len(v) for v in distinct) / len(distinct) > MAX_VALUE_LENGTH:
        return f"values longer than {MAX_VALUE_LENGTH} characters on average (free text)"
    return None


def choose_colour_column(metadata: Metadata, requested: str | None) -> tuple[str | None, str | None]:
    """The column that colours the figures and, when none can, why: the requested one ('none' for no colours;
    a column that does not exist is an error), or the first usable column."""
    if requested is not None:
        if requested.lower() == "none":
            return None, None
        column = _find(metadata.columns, requested.lower())
        if column is None:
            raise BaconError(f"--color-by {requested!r}: no such metadata column (columns: "
                             f"{', '.join(metadata.columns)})")
        reason = unusable_reason(metadata.values(column))
        return (None, f"the metadata column {column!r} cannot colour the figures: {reason}") if reason else \
            (column, None)
    for column in metadata.columns:
        if unusable_reason(metadata.values(column)) is None:
            return column, None
    return None, "no metadata column can colour the figures: " + "; ".join(
        f"{c}: {unusable_reason(metadata.values(c))}" for c in metadata.columns)


def sort_key(value: str) -> tuple:
    """Values in numeric order when they are numbers, otherwise alphabetical without regard to case."""
    try:
        return (0, float(value), value)
    except ValueError:
        return (1, value.lower(), value)


def is_numeric(values: list[str]) -> bool:
    """Whether every present value is a number (and there is one)."""
    present = [v for v in values if v]
    if not present:
        return False
    try:
        return all(f == f for f in map(float, present))  # not NaN
    except ValueError:
        return False


def write_copy(path: Path, metadata: Metadata) -> None:
    """The normalised TSV (sample first, blank cells for missing values), rewritten only when it changes."""
    lines = ["\t".join(["sample", *metadata.columns])]
    lines += ["\t".join([name, *(metadata.value(name, c) for c in metadata.columns)]) for name in metadata.rows]
    text = "\n".join(lines) + "\n"
    if path.exists() and path.read_text(encoding="utf-8") == text:
        return
    tmp = path.with_suffix(".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)
