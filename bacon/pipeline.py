"""The BACoN pipeline: bait, filter, assemble, compare; with checkpoints to resume."""

from __future__ import annotations

import errno
import fcntl
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
            data = _moved(data, root, self.root)
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
    for path in filter(None, outputs):
        for name in FOLDERS.values():
            if f"{os.sep}{name}{os.sep}" in path:
                return path.rsplit(f"{os.sep}{name}{os.sep}", 1)[0]
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


def _output_signature(path: Path | None) -> list[object] | None:
    """Size and modification time of a step's output, recorded by the next step: if the output changes (a sample
    run again, then interrupted), the next step's result is not reused."""
    if path is None or not path.exists():
        return None
    st = path.stat()
    return [st.st_size, st.st_mtime_ns]


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
            return st.sample.name, {"failed": message}
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


def _prepare_reference(s: Settings) -> tuple[Path, int]:
    lengths = check_reference(s.reference)
    s.output.mkdir(parents=True, exist_ok=True)
    local = s.output / "reference.fasta"
    # An uncompressed copy with plain names, upper case and N for ambiguity codes: every tool reads it (ska map is
    # case-sensitive, and SKA2 reads IUPAC codes as fixed bases), and it records what was used.
    records = [Record(r.header, acgtn(r.seq)) for r in read_records(s.reference)]
    tmp = local.with_suffix(".tmp")
    write_fasta(tmp, records)
    if not local.exists() or local.read_bytes() != tmp.read_bytes():
        tmp.replace(local)
    else:
        tmp.unlink()
    return local, sum(n for _, n in lengths)


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
    s.input = s.input.resolve() if s.input else None
    s.sample_sheet = s.sample_sheet.resolve() if s.sample_sheet else None
    s.add_genomes = [p.resolve() for p in s.add_genomes]
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
    if logger.getEffectiveLevel() > logging.INFO:  # The file log always records progress
        logger.setLevel(logging.INFO)
    try:
        return _run(s, started)
    finally:
        logging.getLogger("bacon").removeHandler(file_handler)
        file_handler.close()
        lock.close()  # Releases the lock


def _run(s: Settings, started: float) -> int:
    log.info("BACoN %s", __version__)
    tools = require(needed_tools(s))
    samples = _load_samples(s)
    reference, reference_length = _prepare_reference(s)
    genome_size = s.genome_size or reference_length
    log.info("Reference %s: %s bp in %d sequence(s)%s", s.reference.name, f"{reference_length:,}",
             sum(1 for _ in read_records(reference)),
             f"; genome size set to {genome_size:,} bp" if s.genome_size else "")
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

    def bait(st: SampleState, threads: int) -> steps.StepResult:
        if s.baiting == "minimap2":
            return steps.bait_minimap2(st.sample, reference, folder["bait"], logs / "1_bait", threads, s.keep_bam)
        return steps.bait_bbduk(st.sample, reference, folder["bait"], logs / "1_bait", threads, s.kmer, s.hdist,
                                max(1, s.memory_gb // max(1, min(s.parallel, len(states)))))

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
        # What each sample's step reads: its input files (bait), or the output of the previous step.
        upstream = {st.sample.name: inputs[st.sample.name] if step == "bait" else _output_signature(st.reads)
                    for st in states}
        # Reuse the samples that succeeded with the same parameters and the same input; run the others (new,
        # failed before, or whose input changed). Results of BACoN < 0.3.4 have no "upstream".
        results = {name: res for name, res in (saved or {}).get("results", {}).items()
                   if name in inputs and not res.get("failed") and name not in refreshed
                   and (step != "bait" or res.get("input") == inputs[name])
                   and res.get("upstream", upstream.get(name)) == upstream.get(name)
                   and _outputs_exist({name: res})}
        todo = [st for st in states if not st.failed and st.sample.name not in results]
        if not todo:
            log.info("%s: already done, skipping", labels[step])
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

            new = _run_parallel(todo, functions[step], s, step, logs / LOG_FOLDERS[step], on_done=done)
            refreshed |= {name for name, res in new.items() if not res.get("failed")}
        _apply(states, results, step)
        if all(st.failed for st in states):
            _finish(s, states, tools, None, started)
            raise BaconError(f"All samples failed at the {step} step; see the logs in {logs}")

    for st in states:
        if not st.failed:
            _add_notes(st, genome_size)

    try:
        comparison = _compare(s, states, reference, folder["compare"], logs / "4_compare", checkpoints,
                              fingerprint, bool(refreshed) or redo_from <= STEPS.index("compare"))
    except Exception as exc:  # The assemblies are still worth reporting
        message = str(exc).splitlines()[0] if isinstance(exc, BaconError) else f"{type(exc).__name__}: {exc}"
        _finish(s, states, tools, {"failed": message}, started)
        if isinstance(exc, BaconError):
            raise
        raise BaconError(f"The comparison failed unexpectedly: {message}") from exc
    _finish(s, states, tools, comparison, started)
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
    if saved is not None and saved["results"].pop("inputs", signatures) != signatures:
        saved = None  # An assembly changed since (made again, then interrupted before the comparison)
    if saved is not None and Path(saved["results"].get("distances", "")).is_file():
        log.info("Comparison with %s: already done, skipping", s.snp_method)
        result = saved["results"]
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
            started: float) -> None:
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
        "reference": _reference_info(s),
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
        write_multiqc(s.output, rows, distances, tree)
    except Exception as exc:  # noqa: BLE001
        log.warning("Could not write the MultiQC files: %s", exc)
    try:
        log.info("Report: %s", write_report(s.output))
    except Exception as exc:  # noqa: BLE001
        log.warning("Could not write the HTML report: %s", exc)
    failed = sum(1 for st in states if st.failed)
    log.info("Done: %d sample(s) assembled, %d failed. Results in %s", len(states) - failed, failed, s.output)


def _reference_info(s: Settings) -> dict[str, object]:
    """The reference as given: its path, number of sequences, total length, and the MD5 of the file itself (not
    of BACoN's normalized copy, reference.fasta)."""
    local = s.output / "reference.fasta"
    lengths = [len(r.seq) for r in read_records(local)] if local.exists() else []
    return {"file": str(s.reference), "sequences": len(lengths), "length": sum(lengths),
            "md5": hashlib.md5(s.reference.read_bytes()).hexdigest() if s.reference.is_file() else None}


def default_memory_gb() -> int:
    try:
        total = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):  # pragma: no cover
        return 8
    return max(1, int(total * 0.85 / 1e9))
