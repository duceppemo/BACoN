"""The BACoN pipeline: bait, filter, assemble, compare; with checkpoints to resume."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import platform
import time
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from pathlib import Path

from bacon import BaconError, __version__, compare, steps
from bacon.samples import VALID_NAME, Sample, discover, read_sample_sheet
from bacon.seqio import Record, check_reference, read_records, sniff_format, split_extension, write_fasta
from bacon.tools import require, version

log = logging.getLogger(__name__)

STEPS = ("bait", "filter", "assemble", "compare")
FOLDERS = {"bait": "1_extracted", "filter": "2_filtered", "assemble": "3_assembled", "compare": "4_compared"}
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
        self.folder = output / ".checkpoints"

    def path(self, step: str) -> Path:
        return self.folder / f"{step}.json"

    def load(self, step: str, fingerprint: str) -> dict | None:
        try:
            data = json.loads(self.path(step).read_text())
        except (OSError, ValueError):
            return None
        return data if data.get("fingerprint") == fingerprint else None

    def save(self, step: str, fingerprint: str, results: dict) -> None:
        self.folder.mkdir(parents=True, exist_ok=True)
        tmp = self.path(step).with_suffix(".tmp")
        tmp.write_text(json.dumps({"fingerprint": fingerprint, "results": results}, indent=1, default=str))
        tmp.replace(self.path(step))

    def clear(self, step: str) -> None:
        self.path(step).unlink(missing_ok=True)


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
                  s: Settings, step: str) -> dict[str, dict]:
    """Run `fn` on every sample still in the race; failures of one sample do not stop the others."""
    active = [st for st in states if not st.failed]
    workers = max(1, min(s.parallel, len(active)))
    threads = max(1, s.threads // workers)

    def one(st: SampleState) -> tuple[str, dict]:
        start = time.monotonic()
        try:
            res = fn(st, threads)
        except (steps.SampleFailed, BaconError) as exc:
            message = str(exc).splitlines()[0]
            log.warning("  %s: %s", st.sample.name, message)
            return st.sample.name, {"failed": message}
        log.info("  %s done (%.0f s)", st.sample.name, time.monotonic() - start)
        return st.sample.name, {"output": str(res.output) if res.output else None, "stats": res.stats,
                                "notes": res.notes}

    with ThreadPoolExecutor(max_workers=workers) as pool:
        return dict(pool.map(one, active))


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
    # An uncompressed copy with plain names: every tool reads it, and it records what was used.
    records = list(read_records(s.reference))
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
    row["Note"] = "; ".join(dict.fromkeys(st.notes)) or ""
    return row


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "NA")) for c in columns) + "\n")


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
                 "method": s.baiting, "kmer": s.kmer if s.baiting == "bbduk" else None, "keep_bam": s.keep_bam},
        "filter": {"min_length": s.min_read_length, "keep_percent": s.keep_percent,
                   "target_depth": s.target_depth, "genome_size": genome_size},
        "assemble": {"assembler": s.assembler, "read_type": s.read_type, "min_size": s.min_size,
                     "iterations": s.flye_iterations, "genome_size": genome_size,
                     "template_gaps": s.template_gaps if s.assembler == "samtools" else None},
    }

    def bait(st: SampleState, threads: int) -> steps.StepResult:
        if s.baiting == "minimap2":
            return steps.bait_minimap2(st.sample, reference, folder["bait"], logs / "1_bait", threads, s.keep_bam)
        return steps.bait_bbduk(st.sample, reference, folder["bait"], logs / "1_bait", threads, s.kmer,
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
        # Reuse the samples that succeeded with the same parameters; run the others (new or failed before).
        results = {name: res for name, res in (saved or {}).get("results", {}).items()
                   if name in inputs and not res.get("failed") and name not in refreshed
                   and (step != "bait" or res.get("input") == inputs[name]) and _outputs_exist({name: res})}
        todo = [st for st in states if not st.failed and st.sample.name not in results]
        if not todo:
            log.info("%s: already done, skipping", labels[step])
        else:
            if results:
                log.info("%s: resuming, %d sample(s) already done", labels[step], len(results))
            checkpoints.clear(step)
            log.info("%s...", labels[step])
            new = _run_parallel(todo, functions[step], s, step)
            if step == "bait":
                for name, res in new.items():
                    res["input"] = inputs[name]
            refreshed |= {name for name, res in new.items() if not res.get("failed")}
            results.update(new)
            checkpoints.save(step, fingerprint, results)
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
    except BaconError as exc:  # The assemblies are still worth reporting
        _finish(s, states, tools, {"failed": str(exc).splitlines()[0]}, started)
        raise
    _finish(s, states, tools, comparison, started)
    return 0


def _compare(s: Settings, states: list[SampleState], reference: Path, root: Path, log_dir: Path,
             checkpoints: Checkpoints, fingerprint: str, force: bool) -> dict | None:
    if s.snp_method == "none":
        log.info("Comparison skipped (--snp-method none)")
        return None
    assemblies = {st.sample.name: st.assembly for st in states if not st.failed and st.assembly}
    added = _prepare_added_genomes(s, root, set(assemblies))
    assemblies.update(added)
    if len(assemblies) < 3:
        log.warning("Comparison skipped: a tree needs at least three assemblies (%d available)", len(assemblies))
        return {"skipped": f"only {len(assemblies)} assemblies"}
    fingerprint = _fingerprint(fingerprint, {"method": s.snp_method, "tree": s.tree,
                                             "ska_min_freq": s.ska_min_freq, "assemblies": sorted(assemblies),
                                             "added": {k: _file_signature(v) for k, v in added.items()}})
    saved = None if force else checkpoints.load("compare", fingerprint)
    if saved is not None and Path(saved["results"].get("distances", "")).is_file():
        log.info("Comparison with %s: already done, skipping", s.snp_method)
        result = saved["results"]
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
                                                       min_freq=s.ska_min_freq)
        records = list(read_records(snps))
        names, matrix = compare.snp_distances(records)
        compare.write_distances(out / "snp_distances.tsv", names, matrix)
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
                  "distances": str(out / "snp_distances.tsv"), "alignment": str(tree_input),
                  "core_snps": sites}
        checkpoints.save("compare", fingerprint, result)
    log.info("SNP sites: %s; distances: %s; tree: %s", result["core_snps"], result["distances"],
             result["tree"] or "none")
    return result


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
        write_fasta(target, [Record(f"{name}_{r.name}", r.seq.upper()) for r in read_records(path)])
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
        "samples": {st.sample.name: {"files": [str(f) for f in st.sample.files], "status": st.failed or "ok"}
                    for st in states},
        "comparison": comparison,
    }
    tmp = s.output / "run_info.json.tmp"
    tmp.write_text(json.dumps(info, indent=2, default=str) + "\n")
    tmp.replace(s.output / "run_info.json")
    failed = sum(1 for st in states if st.failed)
    log.info("Done: %d sample(s) assembled, %d failed. Results in %s", len(states) - failed, failed, s.output)


def default_memory_gb() -> int:
    try:
        total = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):  # pragma: no cover
        return 8
    return max(1, int(total * 0.85 / 1e9))
