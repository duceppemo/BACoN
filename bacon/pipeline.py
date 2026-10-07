"""The BACoN pipeline: bait, filter, assemble, compare; with checkpoints to resume."""

from __future__ import annotations

import errno
import fcntl
import gzip
import hashlib
import json
import logging
import os
import platform
import re
import threading
import time
import traceback
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from pathlib import Path

from bacon import BaconError, __version__, compare, steps, tools
from bacon import metadata as md
from bacon.annotation import Annotation, annotation_format, find_regions, genbank_fasta_records, load_annotation
from bacon.multiqc import write_multiqc
from bacon.report import write_report
from bacon.samples import VALID_NAME, Sample, discover, read_sample_sheet
from bacon.seqio import (
    Record,
    acgtn,
    check_reference,
    read_records,
    sniff_format,
    split_extension,
    write_fasta,
)
from bacon.tools import require, version

log = logging.getLogger(__name__)

STEPS = ("bait", "filter", "assemble", "compare")
FOLDERS = {"bait": "1_extracted", "filter": "2_filtered", "assemble": "3_assembled", "compare": "4_compared"}
LOG_FOLDERS = {"bait": "1_bait", "filter": "2_filter", "assemble": "3_assemble", "compare": "4_compare"}
SUMMARY_COLUMNS = [
    "Sample", "Status", "Raw_reads", "Raw_bases", "Baited_reads", "Baited_bases", "Baited_pct",
    "Filtered_reads", "Filtered_bases", "Filtered_N50", "Est_depth",
    "Contigs", "Circular_contigs", "Assembly_length", "Largest_contig", "Assembly_N50", "Length_vs_reference",
    "Assembly_depth", "N_bases", "Note",
]
LOW_DEPTH = 20  # Below this estimated depth of filtered reads, a note is added


@dataclass
class Settings:
    reference: Path
    output: Path
    input: Path | None = None
    sample_sheet: Path | None = None
    annotation: Path | None = None  # Report only: not part of any checkpoint
    metadata: Path | None = None  # Report only
    color_by: str | None = None  # Report only: the metadata column that colours the figures ('none': no colours)
    baiting: str = "minimap2"
    kmer: int = 31
    hdist: int = 1
    keep_bam: bool = False
    min_read_length: int = 500
    keep_percent: float = 95.0
    target_depth: int = 100
    assembler: str = "samtools"
    read_type: str = "nano-hq"
    min_size: int | None = None
    flye_iterations: int = 3
    template_gaps: str = "n"
    genome_size: int | None = None
    snp_method: str = "ska"
    ska_min_freq: float = 1.0
    add_genomes: list[Path] = field(default_factory=list)
    tree: str = "fasttree"
    threads: int = 1
    parallel: int = 2
    memory_gb: int = 8
    redo: str | None = None
    command_line: list[str] = field(default_factory=list)


@dataclass
class SampleState:
    sample: Sample
    reads: Path | None = None
    assembly: Path | None = None
    failed: str = ""
    stats: dict[str, object] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------------------------------------------
# Checkpoints
# ---------------------------------------------------------------------------------------------------------------

def _fingerprint(previous: str, params: dict[str, object]) -> str:
    blob = json.dumps({"previous": previous, "params": params}, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


class Checkpoints:
    """One JSON file per step, recording the parameters it ran with (chained with the previous steps') and
    its per-sample results. A step runs again when its parameters or an earlier step's changed."""

    def __init__(self, output: Path):
        self.root = str(output)
        self.folder = output / ".checkpoints"

    def path(self, step: str) -> Path:
        return self.folder / f"{step}.json"

    def load(self, step: str, fingerprint: str) -> dict | None:
        try:
            data = json.loads(self.path(step).read_text())
        except (OSError, ValueError):
            return None
        if data.get("fingerprint") != fingerprint:
            return None
        root = data.get("root") or _guess_root(data.get("results"))  # No root before 0.3.4
        if root and root != self.root:  # A moved or copied output folder: its own files
            data = {**_moved(data, root, self.root), "relocated": True}
        return data

    def save(self, step: str, fingerprint: str, results: dict) -> None:
        self.folder.mkdir(parents=True, exist_ok=True)
        tmp = self.path(step).with_suffix(".tmp")
        tmp.write_text(json.dumps({"fingerprint": fingerprint, "root": self.root, "results": results}, indent=1,
                                  default=str))
        tmp.replace(self.path(step))

    def clear(self, step: str) -> None:
        self.path(step).unlink(missing_ok=True)


def _guess_root(results: dict | None) -> str | None:
    """The output folder of a checkpoint's results, from the path of an output (not of an input, which can be in
    another BACoN folder) and the names of BACoN's folders."""
    results = results or {}
    outputs = [results.get("distances")] + [r.get("output") for r in results.values() if isinstance(r, dict)]
    for path in (p for p in outputs if isinstance(p, str)):  # "distances" can also be a sample's name
        # The last folder name in the path: a parent folder may have the name of one of BACoN's folders
        at = max(path.rfind(f"{os.sep}{name}{os.sep}") for name in FOLDERS.values())
        if at > 0:
            return path[:at]
    return None


def _moved(data: object, old: str, new: str) -> object:
    """`data` with the paths under the folder `old` moved to `new`."""
    if isinstance(data, str):
        return new + data[len(old):] if data == old or data.startswith(old + os.sep) else data
    if isinstance(data, list):
        return [_moved(x, old, new) for x in data]
    if isinstance(data, dict):
        return {k: _moved(v, old, new) for k, v in data.items()}
    return data


def _same_input(res: dict, upstream: object, reads: Path | None, relocated: bool = False) -> bool:
    """Whether a sample's saved result was made from its current input (`upstream`: the signature of the
    previous step's output; in a moved or copied folder, `relocated`, the sizes alone: a copy may not keep the
    files' times). A result of BACoN < 0.3.5 may have no signature: it must then be newer than its input (an
    input made again after it, by a run interrupted before this step, is newer)."""
    if "upstream" in res:
        return res["upstream"] == upstream or (relocated and _size_only(res["upstream"]) == _size_only(upstream))
    if reads is None or not res.get("output"):
        return True
    try:
        return Path(res["output"]).stat().st_mtime_ns >= reads.stat().st_mtime_ns
    except OSError:
        return False


def _output_signature(path: Path | None) -> list[object] | None:
    """Size and modification time of a step's output, recorded by the next step: if the output changes (a sample
    run again, then interrupted), the next step's result is not reused."""
    if path is None or not path.exists():
        return None
    st = path.stat()
    return [st.st_size, st.st_mtime_ns]


def _size_only(signature: object) -> object:
    """The size of an output signature ([size, mtime_ns], or None): all that a copy of the output folder made
    without the files' times (cp -r, scp, an archive) keeps."""
    return signature[0] if isinstance(signature, list) and signature else signature


def _file_signature(path: Path) -> list[object]:
    st = path.stat()
    return [str(path.resolve()), st.st_size, int(st.st_mtime)]


# ---------------------------------------------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------------------------------------------

def needed_tools(s: Settings) -> list[str]:
    tools = ["minimap2"] if s.baiting == "minimap2" else ["bbduk.sh"]
    if s.keep_bam and s.baiting == "minimap2":
        tools.append("samtools")
    tools += ["filtlong", {"flye": "flye", "myloasm": "myloasm", "samtools": "samtools"}[s.assembler]]
    if s.assembler == "samtools":
        tools.append("minimap2")
    tools += {"parsnp": ["parsnp", "harvesttools"], "ska": ["ska"], "none": []}[s.snp_method]
    if s.snp_method != "none":
        tools.append({"fasttree": "FastTree", "iqtree": "iqtree"}[s.tree])
    return tools


def _load_samples(s: Settings) -> list[Sample]:
    if s.sample_sheet is not None:
        return read_sample_sheet(s.sample_sheet)
    assert s.input is not None
    return discover(s.input)


def _run_parallel(states: list[SampleState], fn: Callable[[SampleState, int], steps.StepResult],
                  s: Settings, step: str, log_dir: Path,
                  on_done: Callable[[str, dict], None] | None = None) -> dict[str, dict]:
    """Run `fn` on every sample still in the race; failures of one sample do not stop the others.

    `on_done(name, result)` is called as each sample finishes (to save progress). On an interruption, the
    programs still running are killed and the samples not started are cancelled."""
    active = [st for st in states if not st.failed]
    workers = max(1, min(s.parallel, len(active)))
    threads = max(1, s.threads // workers)

    def one(st: SampleState) -> tuple[str, dict]:
        start = time.monotonic()
        try:
            res = fn(st, threads)
        except (steps.SampleFailed, BaconError) as exc:
            message = str(exc).splitlines()[0]
            if tools.stopping():  # Killed by the interruption: not a failure of the sample (and not saved)
                log.debug("  %s: %s", st.sample.name, message)
                return st.sample.name, {"failed": "interrupted"}
            log.warning("  %s: %s", st.sample.name, message)
            return st.sample.name, {"failed": message, "stats": getattr(exc, "stats", {})}
        except Exception as exc:  # noqa: BLE001 - an unexpected error fails this sample, not the run
            log_file = log_dir / f"{st.sample.name}.log"
            log_file.parent.mkdir(parents=True, exist_ok=True)
            with open(log_file, "a") as fh:  # The traceback, for a bug report
                fh.write("Unexpected error in BACoN:\n" + traceback.format_exc())
            message = f"unexpected error: {type(exc).__name__}: {exc}".splitlines()[0] + f"; see {log_file}"
            log.warning("  %s: %s", st.sample.name, message)
            return st.sample.name, {"failed": message}
        log.info("  %s done (%.0f s)", st.sample.name, time.monotonic() - start)
        return st.sample.name, {"output": str(res.output) if res.output else None, "stats": res.stats,
                                "notes": res.notes}

    results: dict[str, dict] = {}
    tools.allow_programs()
    pool = ThreadPoolExecutor(max_workers=workers)
    try:
        futures = [pool.submit(one, st) for st in active]
        for future in as_completed(futures):
            name, res = future.result()
            results[name] = res
            if on_done is not None:
                on_done(name, res)
    except BaseException:
        pool.shutdown(wait=False, cancel_futures=True)
        killed = tools.kill_running()
        if killed:
            log.warning("Interrupted: %d running program(s) stopped", killed)
        pool.shutdown(wait=True)
        raise
    pool.shutdown(wait=True)
    return {st.sample.name: results[st.sample.name] for st in active}


def _apply(states: list[SampleState], results: dict[str, dict], step: str) -> None:
    for st in states:
        res = results.get(st.sample.name)
        if st.failed or res is None:
            continue
        if res.get("failed"):
            st.failed = f"failed ({step})"
            st.notes.append(res["failed"])
            st.stats.update(res.get("stats", {}))  # E.g. the raw reads of a sample with none baited
            continue
        st.stats.update(res.get("stats", {}))
        st.notes.extend(res.get("notes", []))
        if res.get("output"):
            if step == "assemble":
                st.assembly = Path(res["output"])
            else:
                st.reads = Path(res["output"])


def _outputs_exist(results: dict[str, dict]) -> bool:
    return all(Path(r["output"]).exists() for r in results.values() if r.get("output"))


def _reference_records(reference: Path) -> list[Record]:
    """The sequences of the reference: a fasta file, or the ORIGIN of a GenBank file (named after the record's
    VERSION, as NCBI's fasta of the same record)."""
    if not reference.is_file():
        raise BaconError(f"Reference file not found: {reference}")
    fmt = annotation_format(reference)
    if fmt == "gff3":
        raise BaconError(f"The reference must be a fasta or GenBank file, not GFF3: {reference} (give the GFF3 "
                         "with --annotation)")
    if fmt == "genbank":
        records = genbank_fasta_records(reference)
        names = [r.name for r in records]
        if len(set(names)) != len(names):
            raise BaconError(f"The reference has duplicate sequence names: {reference}")
        return records
    check_reference(reference)
    return list(read_records(reference))


def _prepare_reference(s: Settings) -> tuple[Path, int]:
    # An uncompressed copy with plain names, upper case and N for ambiguity codes: every tool reads it (ska map is
    # case-sensitive, and SKA2 reads IUPAC codes as fixed bases), and it records what was used.
    records = [Record(r.header, acgtn(r.seq)) for r in _reference_records(s.reference)]
    s.output.mkdir(parents=True, exist_ok=True)
    local = s.output / "reference.fasta"
    tmp = local.with_suffix(".tmp")
    write_fasta(tmp, records)
    if not local.exists() or local.read_bytes() != tmp.read_bytes():
        tmp.replace(local)
    else:
        tmp.unlink()
    return local, sum(len(r.seq) for r in records)


ANNOTATION_COPIES = ("annotation.gb", "annotation.gff3")


def _open_bytes(path: Path):
    """The file's bytes, decompressed when it is gzipped (detected from its content)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rb") if magic == b"\x1f\x8b" else open(path, "rb")


def _record_copy(s: Settings, key: str, name: str | None) -> None:
    """Record in .checkpoints/copies.json the copy (of the annotation or of the metadata) this run writes, or
    None when it has none, as soon as it is written: run_info.json records it too, but only at the end of the
    run, and an interrupted first run would leave BACoN's own copies looking hand-made."""
    path = s.output / ".checkpoints" / "copies.json"
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError):
        data = {}
    data = {**data, key: name} if isinstance(data, dict) else {key: name}
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=1))
    tmp.replace(path)


def _recorded_copy(s: Settings, key: str) -> str | None:
    """The name of the copy (of the annotation or of the metadata) that the previous run of this folder recorded
    as its own, in .checkpoints/copies.json or, for a folder of an earlier version, in run_info.json; None when
    there is none: a file BACoN did not write (hand-made, or from another tool) is kept."""
    for file in (s.output / ".checkpoints" / "copies.json", s.output / "run_info.json"):
        try:
            info = json.loads(file.read_text())
        except (OSError, ValueError):
            continue
        if not isinstance(info, dict) or key not in info:
            continue
        recorded = info[key]
        if file.name == "copies.json":
            return recorded if isinstance(recorded, str) else None
        return recorded.get("copy") if isinstance(recorded, dict) else None
    return None


def _remove_copy(s: Settings, path: Path, recorded: str | None, what: str) -> None:
    """Remove the copy of an earlier run when it was BACoN's; otherwise leave it and say so."""
    if not path.exists():
        return
    if path.name == recorded:
        path.unlink()
    else:
        log.info("%s %s was not written by BACoN: kept (the report uses it; delete it if it is stale)", what, path)


def _reference_regions(records: list[Record], annotation: Annotation | None) -> dict[str, object]:
    """For run_info.json, the LSC/IRb/SSC/IRa regions of each reference sequence (RegionBand.record: from the
    annotation, or from the inverted repeat detected in the sequence; an inverted repeat that is not a plastome
    layout is recorded without regions, with a note), or "none"."""
    regions: dict[str, object] = {}
    for rec in records:
        if annotation is not None and rec.name in annotation.sequences:
            band = annotation.sequences[rec.name].band
        else:
            band = find_regions([], len(rec.seq), rec.seq)
        regions[rec.name] = band.record() if band is not None else "none"
        if band is not None and band.regions:
            log.info("Regions of %s: %s, from %s", rec.name, "/".join(r.name for r in band.regions), band.text())
        elif band is not None and band.note:
            log.info("Regions of %s: none (%s)", rec.name, band.note)
    return regions


def _prepare_annotation(s: Settings, reference: Path) -> tuple[dict[str, object] | None, dict[str, object]]:
    """Validate the annotation (--annotation, or the GenBank reference itself) against the reference's sequences
    and copy it, uncompressed, to OUTPUT/annotation.gb or annotation.gff3 for the report; the copy left by an
    earlier run is removed when this run has no annotation (a file BACoN did not write is kept). Also the
    regions of each reference sequence (_reference_regions), annotation or not."""
    source = s.annotation or (s.reference if annotation_format(s.reference) == "genbank" else None)
    recorded = _recorded_copy(s, "annotation")
    records = list(read_records(reference))
    if source is None:
        for name in ANNOTATION_COPIES:
            _remove_copy(s, s.output / name, recorded, "Annotation")
        _record_copy(s, "annotation", None)
        return None, _reference_regions(records, None)
    sequences = [(r.name, len(r.seq)) for r in records]
    annotation = load_annotation(source, sequences, seqs={r.name: r.seq for r in records})
    for warning in annotation.warnings:
        log.warning("%s", warning)
    copy = s.output / ("annotation.gb" if annotation.format == "genbank" else "annotation.gff3")
    for name in ANNOTATION_COPIES:
        if name != copy.name:
            _remove_copy(s, s.output / name, recorded, "Annotation")
    tmp = copy.with_suffix(".tmp")
    with _open_bytes(source) as src, open(tmp, "wb") as dst:  # The bytes as they are (non-ASCII text included)
        for chunk in iter(lambda: src.read(1 << 20), b""):
            dst.write(chunk)
    if copy.exists() and copy.read_bytes() == tmp.read_bytes():
        tmp.unlink()
    else:
        tmp.replace(copy)
    _record_copy(s, "annotation", copy.name)
    log.info("Annotation %s: %d gene(s) on %d of the reference's %d sequence(s)", source.name, annotation.genes,
             len(annotation.sequences), len(sequences))
    return ({"file": str(source), "format": annotation.format, "copy": copy.name, "genes": annotation.genes,
             "sequences": len(annotation.sequences), "regions": annotation.has_regions,
             "transl_tables": annotation.tables, "md5": hashlib.md5(source.read_bytes()).hexdigest()},
            _reference_regions(records, annotation))


def _prepare_metadata(s: Settings, samples: list[Sample]) -> dict[str, object] | None:
    """Read the metadata (--metadata, and the extra columns of the sample sheet, which --metadata overrides column
    by column), match it to the samples and to the added genomes, choose the column that colours the report, and
    write the normalised copy OUTPUT/metadata.tsv for the report; the copy left by an earlier run is removed when
    this run has no metadata (a metadata.tsv BACoN did not write is kept). Fails before any step on an unreadable
    file, a missing 'sample' column or an unknown --color-by column."""
    given = md.read_metadata(s.metadata) if s.metadata else None
    sheet = md.sheet_metadata(s.sample_sheet) if s.sample_sheet else None
    merged = md.merge(given, sheet)
    copy = s.output / md.COPY_NAME
    if merged is None:
        if s.color_by is not None and s.color_by.lower() != "none":
            raise BaconError(f"--color-by {s.color_by!r}: no metadata (give --metadata, or a sample sheet with "
                             "columns besides 'sample' and 'file')")
        _remove_copy(s, copy, _recorded_copy(s, "metadata"), "Metadata")
        _record_copy(s, "metadata", None)
        return None
    names = [x.name for x in samples]
    added = [parts[0] for p in s.add_genomes if (parts := split_extension(p.name))]  # As _prepare_added_genomes
    metadata, unmatched = md.restrict(merged, names + [a for a in added if a not in names])
    for warning in metadata.warnings:
        log.warning("%s", warning)
    if unmatched:
        log.warning("Metadata: %d row(s) match no sample: %s%s", len(unmatched), ", ".join(unmatched[:5]),
                    " …" if len(unmatched) > 5 else "")
    column, why_not = md.choose_colour_column(metadata, s.color_by)  # BaconError on an unknown column
    if why_not:
        log.warning("Metadata: %s", why_not)
    md.write_copy(copy, metadata)
    _record_copy(s, "metadata", copy.name)
    with_row = sum(1 for n in names if n in merged.rows)
    added_with_row = sum(1 for n in added if n in merged.rows)
    log.info("Metadata: %d column(s) (%s); %d of %d samples have a row%s; %s", len(metadata.columns),
             ", ".join(metadata.columns), with_row, len(names),
             f" (and {added_with_row} of {len(added)} added genomes)" if added else "",
             f"colours by {column}" if column else "no colour column")
    return {"file": str(s.metadata) if s.metadata else None, "copy": copy.name, "columns": metadata.columns,
            "sample_sheet_columns": sheet.columns if sheet else [], "matched": with_row,
            "samples_without_row": len(names) - with_row, "unmatched_rows": len(unmatched), "color_by": column,
            "md5": hashlib.md5(s.metadata.read_bytes()).hexdigest() if s.metadata else None}


def _add_notes(st: SampleState, genome_size: int) -> None:
    bases = st.stats.get("Filtered_bases")
    if isinstance(bases, int):
        depth = bases / genome_size
        st.stats["Est_depth"] = f"{depth:.1f}"
        if depth < LOW_DEPTH:
            st.notes.append(f"low depth ({depth:.0f}x)")
    ratio = st.stats.get("Length_vs_reference")
    if ratio not in (None, "NA") and not 0.8 <= float(ratio) <= 1.2:  # type: ignore[arg-type]
        st.notes.append(f"assembly length is {float(ratio):.2f}x the reference")  # type: ignore[arg-type]


def summary_row(st: SampleState) -> dict[str, str]:
    row = {c: "NA" for c in SUMMARY_COLUMNS}
    row["Sample"] = st.sample.name
    row["Status"] = st.failed or "ok"
    for key, value in st.stats.items():
        if key in row:
            row[key] = str(value)
    # Tabs or line breaks in a message (from a program's error) would break the TSV.
    row["Note"] = " ".join("; ".join(dict.fromkeys(st.notes)).split())
    return row


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "NA")) for c in columns) + "\n")
    tmp.replace(path)


def format_table(columns: list[str], rows: list[dict[str, str]]) -> str:
    widths = [max(len(c), *(len(str(r.get(c, ""))) for r in rows)) for c in columns]
    lines = ["  ".join(c.ljust(w) for c, w in zip(columns, widths))]
    lines += ["  ".join(str(r.get(c, "")).ljust(w) for c, w in zip(columns, widths)) for r in rows]
    return "\n".join(line.rstrip() for line in lines)


def run(s: Settings) -> int:
    started = time.time()
    # Absolute paths: some tools run in their own working folder.
    s.output, s.reference = s.output.resolve(), s.reference.resolve()
    s.annotation = s.annotation.resolve() if s.annotation else None
    s.metadata = s.metadata.resolve() if s.metadata else None
    s.input = s.input.resolve() if s.input else None
    s.sample_sheet = s.sample_sheet.resolve() if s.sample_sheet else None
    s.add_genomes = [p.resolve() for p in s.add_genomes]
    if s.input and s.input.is_dir() and (s.output == s.input or s.output.is_relative_to(s.input)):
        # Its files would be taken for a sample's reads on the next run
        raise BaconError(f"The output folder {s.output} cannot be the input folder or inside it")
    s.output.mkdir(parents=True, exist_ok=True)
    tools.allow_programs()  # After an interruption of an earlier run in this process
    lock = open(s.output / ".bacon.lock", "w")  # noqa: SIM115 - held until the end of the run
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as exc:
        if exc.errno not in (errno.EAGAIN, errno.EACCES):  # A file system without locks (some NFS, SMB mounts)
            log.warning("Could not lock %s (%s): make sure no other BACoN run uses it", s.output, exc.strerror)
        else:
            lock.close()
            raise BaconError(f"Another BACoN run is using {s.output}; wait for it to end, or use another output "
                             "folder") from None
    file_handler = logging.FileHandler(s.output / "bacon.log")
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S"))
    logger = logging.getLogger("bacon")
    logger.addHandler(file_handler)
    level = logger.level
    if logger.getEffectiveLevel() > logging.INFO:  # The file log always records progress
        logger.setLevel(logging.INFO)
    try:
        return _run(s, started)
    except BaconError as exc:  # Logged here, while bacon.log still records
        log.error("%s", exc)
        exc.logged = True  # type: ignore[attr-defined]
        raise
    except KeyboardInterrupt as exc:
        log.error("Interrupted")
        exc.logged = True  # type: ignore[attr-defined]
        raise
    finally:
        logger.removeHandler(file_handler)
        logger.setLevel(level)
        file_handler.close()
        lock.close()  # Releases the lock


def _run(s: Settings, started: float) -> int:
    log.info("BACoN %s", __version__)
    tools = require(needed_tools(s))
    samples = _load_samples(s)
    metadata = _prepare_metadata(s, samples)  # Fails fast on a bad file or column; not part of any checkpoint
    reference, reference_length = _prepare_reference(s)
    genome_size = s.genome_size or reference_length
    log.info("Reference %s: %s bp in %d sequence(s)%s", s.reference.name, f"{reference_length:,}",
             sum(1 for _ in read_records(reference)),
             f"; genome size set to {genome_size:,} bp" if s.genome_size else "")
    annotation, regions = _prepare_annotation(s, reference)  # Fails fast on a bad annotation; not in a checkpoint
    log.info("%d sample(s): %s", len(samples), ", ".join(x.name for x in samples))

    states = [SampleState(x) for x in samples]
    checkpoints = Checkpoints(s.output)
    logs = s.output / "logs"
    folder = {step: s.output / name for step, name in FOLDERS.items()}
    redo_from = STEPS.index(s.redo) if s.redo else len(STEPS)

    # Each sample's input files (path, size, time), checked per sample: adding or changing one sample's reads
    # reruns that sample only.
    inputs = {x.name: [_file_signature(f) for f in x.files] for x in samples}
    params = {
        "bait": {"reference": hashlib.md5(reference.read_bytes()).hexdigest(),
                 "method": s.baiting, "kmer": s.kmer if s.baiting == "bbduk" else None, "keep_bam": s.keep_bam,
                 **({"hdist": s.hdist} if s.baiting == "bbduk" else {})},  # Same checkpoints as 0.3.2 for minimap2
        "filter": {"min_length": s.min_read_length, "keep_percent": s.keep_percent,
                   "target_depth": s.target_depth, "genome_size": genome_size},
        "assemble": {"assembler": s.assembler, "read_type": s.read_type, "min_size": s.min_size,
                     "iterations": s.flye_iterations, "genome_size": genome_size,
                     "template_gaps": s.template_gaps if s.assembler == "samtools" else None},
    }

    baiting = {"samples": len(states)}  # Samples to bait in this run

    def bait(st: SampleState, threads: int) -> steps.StepResult:
        if s.baiting == "minimap2":
            return steps.bait_minimap2(st.sample, reference, folder["bait"], logs / "1_bait", threads, s.keep_bam)
        # The memory is shared between the samples baited at the same time
        return steps.bait_bbduk(st.sample, reference, folder["bait"], logs / "1_bait", threads, s.kmer, s.hdist,
                                max(1, s.memory_gb // max(1, min(s.parallel, baiting["samples"]))))

    def filt(st: SampleState, threads: int) -> steps.StepResult:
        assert st.reads is not None
        return steps.filter_filtlong(st.sample.name, st.reads, st.sample.fmt, folder["filter"],
                                     logs / "2_filter", min_length=s.min_read_length,
                                     keep_percent=s.keep_percent, target_bases=genome_size * s.target_depth)

    dirs = steps.AssemblyDirs(folder["assemble"])

    def assemble(st: SampleState, threads: int) -> steps.StepResult:
        assert st.reads is not None
        name, reads, log_dir = st.sample.name, st.reads, logs / "3_assemble"
        common = {"threads": threads, "reference_length": reference_length}
        if s.assembler == "flye":
            return steps.assemble_flye(name, reads, dirs, log_dir, genome_size=genome_size, read_type=s.read_type,
                                       min_overlap=s.min_size, iterations=s.flye_iterations, **common)
        if s.assembler == "myloasm":
            return steps.assemble_myloasm(name, reads, dirs, log_dir, **common)
        return steps.assemble_samtools(name, reads, reference, dirs, log_dir,
                                       fill_gaps=s.template_gaps == "reference", **common)

    functions = {"bait": bait, "filter": filt, "assemble": assemble}
    labels = {"bait": f"Baiting reads matching the reference with {s.baiting}",
              "filter": "Filtering reads with Filtlong",
              "assemble": {"samtools": "Building the templated consensus with samtools",
                           "flye": "Assembling with Flye", "myloasm": "Assembling with myloasm"}[s.assembler]}
    fingerprint = __version__.split(".")[0]
    refreshed: set[str] = set()  # Samples whose output changed in this run: their later steps must run again
    for i, step in enumerate(STEPS[:-1]):
        fingerprint = _fingerprint(fingerprint, params[step])
        saved = None if i >= redo_from else checkpoints.load(step, fingerprint)
        relocated = bool(saved and saved.get("relocated"))  # A moved or copied folder: sizes alone are compared
        # What each sample's step reads: its input files (bait), or the output of the previous step.
        upstream = {st.sample.name: inputs[st.sample.name] if step == "bait" else _output_signature(st.reads)
                    for st in states}
        reads = {st.sample.name: st.reads for st in states}
        # Reuse the samples that succeeded with the same parameters and the same input; run the others (new,
        # failed before, or whose input changed).
        results = {name: res for name, res in (saved or {}).get("results", {}).items()
                   if name in inputs and not res.get("failed") and name not in refreshed
                   and (res.get("input") == inputs[name] if step == "bait"
                        else _same_input(res, upstream.get(name), reads.get(name), relocated))
                   and _outputs_exist({name: res})}
        for name, res in results.items():  # The input's signature, as it is here: a copy's files, earlier versions
            if step != "bait":
                res["upstream"] = upstream[name]
        todo = [st for st in states if not st.failed and st.sample.name not in results]
        if not todo:
            log.info("%s: already done, skipping", labels[step])
            checkpoints.save(step, fingerprint, results)  # With the signatures and the folder of this version
        else:
            if results:
                log.info("%s: resuming, %d sample(s) already done", labels[step], len(results))
            log.info("%s...", labels[step])
            # Saved at once: a checkpoint left by other parameters must not outlive the start of the step (an
            # interruption could leave outputs made with these parameters under it).
            checkpoints.save(step, fingerprint, results)
            lock = threading.Lock()

            def done(name: str, res: dict, step: str = step, fp: str = fingerprint, results: dict = results,
                     upstream: dict = upstream, lock: threading.Lock = lock) -> None:
                # Saved as each sample finishes: an interruption loses only the samples still running.
                if step == "bait":
                    res["input"] = inputs[name]
                elif not res.get("failed"):
                    res["upstream"] = upstream[name]
                with lock:
                    results[name] = res
                    checkpoints.save(step, fp, results)

            if step == "bait":
                baiting["samples"] = len(todo)
            new = _run_parallel(todo, functions[step], s, step, logs / LOG_FOLDERS[step], on_done=done)
            refreshed |= {name for name, res in new.items() if not res.get("failed")}
        _apply(states, results, step)
        if all(st.failed for st in states):
            _finish(s, states, tools, None, started, annotation, metadata, regions)
            raise BaconError(f"All samples failed at the {step} step; see the logs in {logs}")

    for st in states:
        if not st.failed:
            _add_notes(st, genome_size)

    try:
        comparison = _compare(s, states, reference, folder["compare"], logs / "4_compare", checkpoints,
                              fingerprint, bool(refreshed) or redo_from <= STEPS.index("compare"))
    except Exception as exc:  # The assemblies are still worth reporting
        message = str(exc).splitlines()[0] if isinstance(exc, BaconError) else f"{type(exc).__name__}: {exc}"
        _finish(s, states, tools, {"failed": message}, started, annotation, metadata, regions)
        if isinstance(exc, BaconError):
            raise
        raise BaconError(f"The comparison failed unexpectedly: {message}") from exc
    _finish(s, states, tools, comparison, started, annotation, metadata, regions)
    return 0


def _compare(s: Settings, states: list[SampleState], reference: Path, root: Path, log_dir: Path,
             checkpoints: Checkpoints, fingerprint: str, force: bool) -> dict | None:
    if s.snp_method == "none":
        log.info("Comparison skipped (--snp-method none)")
        return None
    assemblies = {st.sample.name: st.assembly for st in states if not st.failed and st.assembly}
    added = _prepare_added_genomes(s, root, {st.sample.name for st in states})  # Failed samples' names too
    assemblies.update(added)
    if len(assemblies) < 3:
        log.warning("Comparison skipped: a tree needs at least three assemblies (%d available)", len(assemblies))
        return {"skipped": f"only {len(assemblies)} assemblies"}
    fingerprint = _fingerprint(fingerprint, {"method": s.snp_method, "tree": s.tree,
                                             "ska_min_freq": s.ska_min_freq, "assemblies": sorted(assemblies),
                                             "added": {k: hashlib.md5(v.read_bytes()).hexdigest()
                                                       for k, v in added.items()}})
    signatures = {k: _output_signature(v) for k, v in sorted(assemblies.items())}
    saved = None if force else checkpoints.load("compare", fingerprint)
    if saved is not None and not _same_assemblies(saved["results"].pop("inputs", None), signatures, assemblies,
                                                  saved["results"].get("distances"), bool(saved.get("relocated"))):
        saved = None  # An assembly changed since (made again, then interrupted before the comparison)
    if saved is not None and Path(saved["results"].get("distances", "")).is_file():
        log.info("Comparison with %s: already done, skipping", s.snp_method)
        result = saved["results"]
        checkpoints.save("compare", fingerprint, {**result, "inputs": signatures})  # Signatures of this version
        if not result.get("vcf") or not Path(result["vcf"]).is_file() or _vcf_outdated(Path(result["vcf"])):
            # A comparison from BACoN < 0.3.2, whose VCF failed, or whose VCF predates the fixes of 0.3.3: the VCF
            # is made from the comparison's files.
            out = Path(result["distances"]).parent
            paths = {k: v for k, v in assemblies.items() if v is not None}
            result["vcf"] = _write_vcf(s, reference, out, log_dir, paths, set(added))
            checkpoints.save("compare", fingerprint, {**result, "inputs": signatures})
    else:
        checkpoints.clear("compare")
        out = root / (s.snp_method if s.snp_method != "ska" or s.ska_min_freq == 1 else f"ska_{s.ska_min_freq:g}")
        log.info("Comparing %d assemblies with %s...", len(assemblies),
                 "SKA2" if s.snp_method == "ska" else "Parsnp")
        paths = {k: v for k, v in assemblies.items() if v is not None}
        if s.snp_method == "parsnp":
            tree_input, snps = compare.run_parsnp(reference, paths, out, log_dir, threads=s.threads)
        else:
            tree_input, snps = compare.run_ska(reference, paths, out, log_dir, threads=s.threads,
                                               min_freq=s.ska_min_freq, added=set(added))
        records = list(read_records(snps))
        names, matrix = compare.snp_distances(records)
        compare.write_distances(out / "snp_distances.tsv", names, matrix)
        vcf = _write_vcf(s, reference, out, log_dir, paths, set(added))
        sites = _alignment_length(snps)
        tree = None
        if sites == 0:
            hint = (" Some assemblies may be incomplete; try --ska-min-freq below 1 (e.g. 0.9) to keep SNPs "
                    "missing from a few genomes." if s.snp_method == "ska" else "")
            log.warning("No SNP site is shared by all the genomes: no tree.%s", hint)
        else:
            log.info("Building a tree with %s...", "FastTree" if s.tree == "fasttree" else "IQ-TREE")
            tree = str(compare.build_tree(tree_input, out, log_dir, method=s.tree, threads=s.threads))
        result = {"method": s.snp_method, "tree_method": s.tree, "tree": tree,
                  "tree_svg": str(out / "tree.svg") if tree else None,
                  "distances": str(out / "snp_distances.tsv"), "alignment": str(tree_input), "vcf": vcf,
                  "core_snps": sites}
        checkpoints.save("compare", fingerprint, {**result, "inputs": signatures})
    log.info("SNP sites: %s; distances: %s; tree: %s", result["core_snps"], result["distances"],
             result["tree"] or "none")
    return result


def _same_assemblies(saved: dict | None, signatures: dict, assemblies: dict[str, Path | None],
                     distances: str | None, relocated: bool = False) -> bool:
    """Whether a saved comparison was made from the current assemblies (in a moved or copied folder,
    `relocated`, by their sizes alone). A comparison of BACoN < 0.3.4 has no signatures: it must then be newer
    than every assembly."""
    if saved is not None:
        return saved == signatures or (relocated and {k: _size_only(v) for k, v in saved.items()}
                                       == {k: _size_only(v) for k, v in signatures.items()})
    try:
        made = Path(distances or "").stat().st_mtime_ns
        return all(p is None or p.stat().st_mtime_ns <= made for p in assemblies.values())
    except OSError:
        return False


# VCFs written by earlier versions are rewritten on resume: 0.3.3 fixed N alleles, SNPs next to the ends of
# circular genomes and Parsnp's INFO values.
VCF_REWRITE_BEFORE = (0, 3, 3)


def _vcf_outdated(vcf: Path) -> bool:
    """True if the VCF was written by a BACoN older than VCF_REWRITE_BEFORE (or does not say)."""
    with open(vcf, encoding="ascii", errors="replace") as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            if line.startswith("##source=BACoN "):
                numbers = [int(n) for n in re.findall(r"\d+", line[len("##source=BACoN "):])[:3]]
                return tuple(numbers) < VCF_REWRITE_BEFORE  # () if no version
    return True


def _write_vcf(s: Settings, reference: Path, out: Path, log_dir: Path, paths: dict[str, Path],
               added: set[str]) -> str | None:
    """snps.vcf of the comparison; an extra output, whose failure is logged and does not fail the comparison."""
    try:
        vcf_path, records = compare.write_vcf(s.snp_method, reference, out, log_dir, threads=s.threads,
                                              assemblies=paths, source=f"BACoN {__version__}", added=added)
    except Exception as exc:  # noqa: BLE001
        log.warning("Could not write the VCF: %s", str(exc).splitlines()[0] if str(exc) else repr(exc))
        return None
    log.info("VCF: %s (%d SNP records)", vcf_path, records)
    return str(vcf_path)


def _prepare_added_genomes(s: Settings, root: Path, samples: set[str]) -> dict[str, Path]:
    """Genomes given with --add-genomes, copied as plain fasta for the comparison, named after their files."""
    added: dict[str, Path] = {}
    for path in s.add_genomes:
        parts = split_extension(path.name)
        if not parts or not path.is_file() or sniff_format(path) != "fasta":
            raise BaconError(f"--add-genomes: not a fasta file: {path}")
        name = parts[0]
        if not VALID_NAME.match(name) or name.lower() == "reference" or name in samples or name in added:
            raise BaconError(f"--add-genomes: the name {name!r} (from {path}) is invalid or already used; "
                             "rename the file")
        target = root / "added_genomes" / f"{name}.fasta"
        target.parent.mkdir(parents=True, exist_ok=True)
        tmp = target.with_suffix(".tmp")
        write_fasta(tmp, [Record(f"{name}_{r.name}", acgtn(r.seq)) for r in read_records(path)])
        if target.exists() and target.read_bytes() == tmp.read_bytes():
            tmp.unlink()  # Unchanged: the copy (and the comparison's fingerprint) stays as it was
        else:
            tmp.replace(target)
        added[name] = target
    if added:
        log.info("Added to the comparison: %s", ", ".join(added))
    return added


def _alignment_length(path: Path) -> int:
    for rec in read_records(path):
        return len(rec.seq)
    return 0


def _finish(s: Settings, states: list[SampleState], tools: dict[str, str], comparison: dict | None,
            started: float, annotation: dict[str, object] | None = None,
            metadata: dict[str, object] | None = None, regions: dict[str, object] | None = None) -> None:
    rows = [summary_row(st) for st in states]
    write_tsv(s.output / "summary.tsv", SUMMARY_COLUMNS, rows)
    shown = ["Sample", "Status", "Baited_reads", "Filtered_reads", "Est_depth", "Contigs", "Assembly_length",
             "Note"]
    log.info("Summary (%s):\n%s", s.output / "summary.tsv", format_table(shown, rows))
    info = {
        "bacon_version": __version__,
        "command_line": s.command_line,
        "started": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(started)),
        "duration_s": round(time.time() - started, 1),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "settings": {k: (str(v) if isinstance(v, Path) else v) for k, v in asdict(s).items()
                     if k != "command_line"},
        "tools": {name: {"path": path, "version": version(name)} for name, path in tools.items()},
        "reference": _reference_info(s, regions),
        "annotation": annotation,
        "metadata": metadata,
        "samples": {st.sample.name: {"files": [str(f) for f in st.sample.files], "status": st.failed or "ok"}
                    for st in states},
        "comparison": comparison,
    }
    tmp = s.output / "run_info.json.tmp"
    tmp.write_text(json.dumps(info, indent=2, default=str) + "\n")
    tmp.replace(s.output / "run_info.json")
    # The report and the MultiQC files are conveniences: their failure must not fail the run.
    try:
        distances = Path(comparison["distances"]) if comparison and comparison.get("distances") else None
        tree = Path(comparison["tree"]) if comparison and comparison.get("tree") else None
        columns = md.read_metadata(s.output / md.COPY_NAME) if metadata else None
        write_multiqc(s.output, rows, distances, tree, columns)
    except Exception as exc:  # noqa: BLE001
        log.warning("Could not write the MultiQC files: %s", exc)
    try:
        log.info("Report: %s", write_report(s.output))
    except Exception as exc:  # noqa: BLE001
        log.warning("Could not write the HTML report: %s", exc)
    failed = sum(1 for st in states if st.failed)
    log.info("Done: %d sample(s) assembled, %d failed. Results in %s", len(states) - failed, failed, s.output)


def _reference_info(s: Settings, regions: dict[str, object] | None = None) -> dict[str, object]:
    """The reference as given: its path, number of sequences, total length, the MD5 of the file itself (fasta or
    GenBank; not of BACoN's normalized copy, reference.fasta), and the regions of each sequence (or "none")."""
    local = s.output / "reference.fasta"
    lengths = [len(r.seq) for r in read_records(local)] if local.exists() else []
    return {"file": str(s.reference), "sequences": len(lengths), "length": sum(lengths),
            "md5": hashlib.md5(s.reference.read_bytes()).hexdigest() if s.reference.is_file() else None,
            "regions": regions or {}}


CGROUP_MEMORY_FILES = ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes")  # v2, v1


def _cgroup_memory_limit(files: tuple[str, ...] = CGROUP_MEMORY_FILES) -> int | None:
    """The memory limit of this process's control group (a job scheduler's or a container's), in bytes."""
    for path in files:
        try:
            value = Path(path).read_text().strip()
        except OSError:
            continue
        if value.isdigit() and int(value) < 1 << 60:  # "max" or a huge number: no limit
            return int(value)
    return None


def default_memory_gb() -> int:
    """85% of the memory BACoN may use: the physical memory, or less under a job scheduler or in a container."""
    try:
        total = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):  # pragma: no cover
        return 8
    total = min(total, _cgroup_memory_limit() or total)
    return max(1, int(total * 0.85 / 1e9))


def usable_cpus() -> int:
    """The CPUs BACoN may use: all of them, or those a job scheduler or `taskset` gave it."""
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0)) or 1
    return os.cpu_count() or 1  # pragma: no cover - macOS
