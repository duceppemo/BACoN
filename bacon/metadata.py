"""Sample metadata (`--metadata`, or the extra columns of a sample sheet): read, merged, and the column that
colours the report's figures. Metadata is for the report only: it is not part of any checkpoint."""

from __future__ import annotations

import csv
import logging
import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

from bacon import BaconError

log = logging.getLogger(__name__)

MISSING_VALUES = {"", "NA", "na", "-"}  # Cells read as "no value"
PALETTE_COLOURS = 12  # The report's categorical colours (bacon.report.PALETTE_LIGHT and PALETTE_DARK)
# Distinct values a column may have to colour the figures: each value has a marker of its own, a colour and a shape
# (bacon.report.MARKERS)
MAX_COLOUR_VALUES = 48
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


def read_table(path: Path, what: str, comment_lines: bool = False,
               encoding_warning: bool = True) -> tuple[list[str], list[dict[str, str]]]:
    """Header and rows of a TSV or CSV file: UTF-8 with or without BOM, UTF-16 with a BOM, or else Windows-1252
    (with a warning, unless not `encoding_warning`); LF, CRLF or CR line ends (only those: other characters that
    Python counts as line breaks, such as a form feed, are whitespace in a cell). Blank lines are ignored, and
    lines starting with '#' before the header are comments (after it, they are data, unless `comment_lines`: a
    sample sheet). A tab in the header makes it a TSV: cells are split on tabs only (no quote can swallow a row),
    and a cell entirely in quotes loses them (`""` inside is one quote); otherwise a CSV, whose quoted values may
    hold commas but not line breaks. Names and values are stripped (names also lose their inner whitespace);
    columns with an empty name (a trailing separator) are dropped, and two columns of the same name (any case)
    are an error."""
    if not path.is_file():
        raise BaconError(f"{what} not found: {path}")
    try:
        data = path.read_bytes()
    except OSError as exc:
        raise BaconError(f"{what} cannot be read: {path} ({exc.strerror})") from None
    text = _decode(data, path, what, encoding_warning)
    lines: list[tuple[int, str]] = []  # (line number in the file, text)
    for number, line in enumerate(text.replace("\r\n", "\n").replace("\r", "\n").split("\n"), 1):
        if not line.strip() or (line.lstrip().startswith("#") and (comment_lines or not lines)):
            continue
        lines.append((number, line))
    if not lines:
        raise BaconError(f"{what} is empty: {path}")
    tsv = "\t" in lines[0][1]
    records = ([[_unquoted(c) for c in line.split("\t")] for _, line in lines] if tsv
               else _csv_records(lines, path, what))
    header = [_squash(c) for c in records[0]]
    lower = [c.lower() for c in header if c]
    duplicates = [c for c in dict.fromkeys(header) if c and lower.count(c.lower()) > 1]  # In header order
    if duplicates:
        raise BaconError(f"{what} {path}: the column name {', '.join(repr(d) for d in duplicates)} is used more "
                         "than once (names differing only in case are the same): give each column its own name")
    rows = []
    for cells in records[1:]:
        cells += [""] * (len(header) - len(cells))  # Fewer cells than names: empty; more: ignored
        rows.append({name: value.strip() for name, value in zip(header, cells) if name})
    return [c for c in header if c], rows


def _decode(data: bytes, path: Path, what: str, warn: bool) -> str:
    """The text of a table: UTF-16 when it starts with a UTF-16 byte order mark (as some spreadsheets save
    "Unicode text"), UTF-8 (with or without a BOM), or Windows-1252 (Excel's encoding on Windows) when it is not
    valid UTF-8."""
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        return data.decode("utf-16", errors="replace")
    try:
        return data.decode("utf-8-sig")
    except UnicodeDecodeError:
        if warn:
            log.warning("%s %s is not UTF-8: read as Windows-1252 (cp1252); if some characters look wrong, save "
                        "it as UTF-8", what, path)
        return data.decode("cp1252", errors="replace")


_END = "\n"  # After the last line: a quoted value that reaches it was never closed


def _csv_records(lines: list[tuple[int, str]], path: Path, what: str) -> list[list[str]]:
    """The cells of each line of a CSV, as 0.3.5 read them (Python's csv module, not strict: `"a" ,b` is `a `
    and b), but an unclosed quote or a quoted line break is an error naming the line (the value would otherwise
    swallow the following rows, or break the copy)."""
    reader = csv.reader([*(line for _, line in lines), _END])
    records = []
    while True:
        start = reader.line_num
        try:
            cells = next(reader)
        except StopIteration:
            return records
        except csv.Error as exc:
            raise BaconError(f"{what} {path}, line {lines[start][0]}: {exc}") from None
        if reader.line_num > len(lines):  # Into the end marker: the quote was never closed
            if start == len(lines):  # The end marker itself
                return records
            raise BaconError(f"{what} {path}, line {lines[start][0]}: an unclosed quote (its value would swallow "
                             "the lines after it)")
        if reader.line_num - start > 1:
            raise BaconError(f"{what} {path}, line {lines[start][0]}: a quoted value spans several lines (an "
                             "unclosed quote?)")
        records.append(cells)


def _unquoted(cell: str) -> str:
    """A TSV cell without the pair of quotes around it, when it is entirely quoted (as a spreadsheet exports
    it, and as 0.3.5 read it): `"abc"` is abc, `""` empty, `"a""b"` a"b; `5" tube` is left as it is."""
    cell = cell.strip()
    if len(cell) >= 2 and cell[0] == cell[-1] == '"':
        return cell[1:-1].replace('""', '"')
    return cell


def _squash(value: str) -> str:
    """Stripped, inner runs of whitespace (tabs, line breaks) as one space: safe in the TSV copy."""
    return " ".join(value.split())


def _clean(value: str) -> str:
    """A metadata value: missing values as "", the others without stray whitespace."""
    return "" if value in MISSING_VALUES else _squash(value)


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
    # Its encoding was reported when its samples were read (bacon.samples.read_sample_sheet)
    header, rows = read_table(path, "Sample sheet", comment_lines=True, encoding_warning=False)
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
    lacks; names differing only in case are the same column (under the first's name)."""
    if first is None or second is None:
        return first or second
    taken = {c.lower() for c in first.columns}
    columns = first.columns + [c for c in second.columns if c.lower() not in taken]
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
        return (f"{len(distinct)} distinct values (at most {MAX_COLOUR_VALUES} can be coloured, each with a "
                "marker of its own)")
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


NUMBER = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$", re.ASCII)  # What the page's parseFloat reads whole
_CHUNKS = re.compile(r"(\d+)", re.ASCII)


def sort_key(value: str) -> tuple:
    """Values in numeric order when they are numbers, otherwise in natural order (v2 before v10) without regard to
    case, as the report's table sorts them."""
    if NUMBER.match(value):
        return (0, float(value), value)
    natural = tuple((0, int(p), "") if p.isdigit() else (1, 0, p.lower()) for p in _CHUNKS.split(value) if p)
    return (1, natural, value)


def is_numeric(values: list[str]) -> bool:
    """Whether every present value is a plain number (digits, a decimal point, an exponent; and there is one)."""
    present = [v for v in values if v]
    return bool(present) and all(NUMBER.match(v) for v in present)


def shown_name(column: str, taken: Iterable[str]) -> str:
    """The column's name in a table that has columns of its own: 'NAME (metadata)' when it is one of theirs."""
    return f"{column} (metadata)" if column.lower() in {t.lower() for t in taken} else column


def _tsv_cell(value: str) -> str:
    """A value as the copy holds it: a value in quotes ("A") is quoted again, so that read_table, which removes the
    quotes around a cell, gives it back."""
    if len(value) >= 2 and value[0] == value[-1] == '"':
        return '"' + value.replace('"', '""') + '"'
    return value


def write_copy(path: Path, metadata: Metadata) -> None:
    """The normalised TSV (sample first, blank cells for missing values), rewritten only when it changes; read
    with read_metadata, it gives the same values."""
    rows = [["sample", *metadata.columns]]
    rows += [[name, *(metadata.value(name, c) for c in metadata.columns)] for name in metadata.rows]
    data = ("\n".join("\t".join(_tsv_cell(v) for v in row) for row in rows) + "\n").encode("utf-8")
    try:
        if path.read_bytes() == data:  # Bytes: a file in another encoding (hand-made) is simply replaced
            return
    except OSError:
        pass
    tmp = path.with_suffix(".tmp")
    tmp.write_bytes(data)
    tmp.replace(path)
