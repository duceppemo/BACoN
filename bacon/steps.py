"""Per-sample steps: bait, filter and assemble reads."""

from __future__ import annotations

import logging
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from bacon import BaconError
from bacon.samples import Sample
from bacon.seqio import (
    ReadStats,
    Record,
    concatenate,
    extract_reads,
    n50,
    open_write,
    read_records,
    read_stats,
    write_fasta,
)
from bacon.tools import run, which

log = logging.getLogger(__name__)


class SampleFailed(BaconError):
    """This sample cannot go further; the others can."""


@dataclass
class StepResult:
    output: Path | None
    stats: dict[str, object] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


def reads_suffix(fmt: str) -> str:
    return ".fastq.gz" if fmt == "fastq" else ".fasta.gz"


# ---------------------------------------------------------------------------------------------------------------
# 1. Baiting
# ---------------------------------------------------------------------------------------------------------------

def _names_from_paf(paf: Path) -> set[str]:
    with open(paf) as fh:
        return {line.split("\t", 1)[0] for line in fh if line.strip()}


def bait_minimap2(sample: Sample, reference: Path, out_dir: Path, log_dir: Path, threads: int,
                  keep_bam: bool) -> StepResult:
    """Keep the reads with at least one alignment to the reference."""
    out_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{sample.name}.log"
    mm2 = ["minimap2", "-x", "map-ont", "--secondary=no", "-t", str(threads)]
    files = [str(f) for f in sample.files]
    if keep_bam:
        bam = out_dir / f"{sample.name}.bam"
        run([[*mm2, "-a", str(reference), *files],
             ["samtools", "sort", "-@", str(max(1, threads // 4)), "-o", str(bam), "-"]], log_file,
            what=f"(mapping {sample.name})")
        run(["samtools", "index", str(bam)], log_file)
        names_file = out_dir / f".{sample.name}.names"
        run(["samtools", "view", "-F", "0x4", str(bam)], log_file, stdout=names_file)
        with open(names_file) as fh:
            names = {line.split("\t", 1)[0] for line in fh if line.strip()}
        names_file.unlink()
    else:
        paf = out_dir / f".{sample.name}.paf"
        run([*mm2, str(reference), *files], log_file, stdout=paf, what=f"(mapping {sample.name})")
        names = _names_from_paf(paf)
        paf.unlink()
    output = out_dir / (sample.name + reads_suffix(sample.fmt))
    total, kept = extract_reads(sample.files, names, output)
    return _bait_result(output, total, kept)


_BBDUK_COUNTS = re.compile(r"^(Input|Contaminants):\s+(\d+) reads\s.*?(\d+) bases", re.M)


def bait_bbduk(sample: Sample, reference: Path, out_dir: Path, log_dir: Path, threads: int, kmer: int,
               memory_gb: int) -> StepResult:
    """Keep the reads sharing at least one k-mer (up to 2 mismatches) with the reference."""
    out_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{sample.name}.log"
    source = sample.files[0]
    merged = None
    if len(sample.files) > 1:  # BBDuk takes one input file
        merged = out_dir / f".{sample.name}.merged{reads_suffix(sample.fmt)}"
        concatenate(sample.files, merged)
        source = merged
    output = out_dir / (sample.name + reads_suffix(sample.fmt))
    start = log_file.stat().st_size if log_file.exists() else 0
    try:
        run(["bbduk.sh", f"-Xmx{max(1, memory_gb)}g", "-eoom", "overwrite=true", f"in={source}",
             f"ref={reference}", f"threads={threads}", f"k={kmer}", "hdist=2", "maskmiddle=f",
             f"outm={output}"], log_file, what=f"(baiting {sample.name})")
    finally:
        if merged:
            merged.unlink(missing_ok=True)
    with open(log_file) as fh:
        fh.seek(start)
        counts = {m.group(1): (int(m.group(2)), int(m.group(3))) for m in _BBDUK_COUNTS.finditer(fh.read())}
    total = ReadStats(*counts["Input"]) if "Input" in counts else None
    return _bait_result(output, total, read_stats([output]))


def _bait_result(output: Path, total: ReadStats | None, kept: ReadStats) -> StepResult:
    stats: dict[str, object] = {"Raw_reads": total.reads if total else "NA",
                                "Raw_bases": total.bases if total else "NA",
                                "Baited_reads": kept.reads, "Baited_bases": kept.bases}
    stats["Baited_pct"] = (f"{100 * kept.bases / total.bases:.3f}" if total and total.bases else "NA")
    if kept.reads == 0:
        raise SampleFailed("no reads matched the reference")
    return StepResult(output, stats)


# ---------------------------------------------------------------------------------------------------------------
# 2. Read filtering
# ---------------------------------------------------------------------------------------------------------------

def filter_filtlong(name: str, reads: Path, fmt: str, out_dir: Path, log_dir: Path, *, min_length: int,
                    keep_percent: float, target_bases: int) -> StepResult:
    """Filtlong for fastq reads; for fasta reads (no qualities, which Filtlong needs), `filter_by_length`."""
    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / (name + reads_suffix(fmt))
    if fmt == "fasta":
        filter_by_length(reads, output, min_length=min_length, target_bases=target_bases)
    else:
        run(["filtlong", "--min_length", str(min_length), "--keep_percent", f"{keep_percent:g}",
             "--target_bases", str(target_bases), str(reads)], log_dir / f"{name}.log", stdout=output,
            what=f"(filtering {name})")
    stats = read_stats([output])
    if stats.reads == 0:
        raise SampleFailed(f"no reads left after filtering (all shorter than {min_length} bp?)")
    return StepResult(output, {"Filtered_reads": stats.reads, "Filtered_bases": stats.bases,
                               "Filtered_N50": stats.n50})


def filter_by_length(reads: Path, output: Path, *, min_length: int, target_bases: int) -> None:
    """Keep the longest reads of at least `min_length` bp, up to `target_bases` (in their original order)."""
    lengths = sorted((len(r.seq) for r in read_records(reads) if len(r.seq) >= min_length), reverse=True)
    total, threshold, ties = 0, 0, 0
    for length in lengths:  # The shortest length still kept, and how many reads of that length to keep
        if total >= target_bases:
            break
        total += length
        ties = ties + 1 if length == threshold else 1
        threshold = length
    with open_write(output) as out:
        for rec in read_records(reads):
            n = len(rec.seq)
            if n > threshold or (n == threshold and n >= min_length and ties > 0):
                if n == threshold:
                    ties -= 1
                out.write(rec.format())


# ---------------------------------------------------------------------------------------------------------------
# 3. Assembly
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class AssemblyDirs:
    root: Path  # 3_assembled

    @property
    def assemblies(self) -> Path:
        return self.root / "all_assemblies"

    @property
    def graphs(self) -> Path:
        return self.root / "assembly_graphs"

    def work(self, name: str) -> Path:
        return self.root / name


def _publish_assembly(name: str, contigs: Path, dirs: AssemblyDirs, *, rename: bool,
                      circular: set[str] | None = None) -> list[Record]:
    """Copy the contigs to all_assemblies/<sample>.fasta, prefixing the contig names with the sample name.

    Contigs known to be circular (in `circular`, or flagged so in their header by the assembler) are marked
    `circular=true` in their new header.
    """
    records = [r for r in read_records(contigs) if r.seq]
    if not records:
        raise SampleFailed("the assembler produced no contigs")
    dirs.assemblies.mkdir(parents=True, exist_ok=True)
    out = []
    for i, rec in enumerate(records, 1):
        header = f"{name}_{i}" if rename else f"{name}_{rec.name}"
        if rec.name in (circular or set()) or _CIRCULAR_TAGS.search(rec.header):
            header += " circular=true"
        out.append(Record(header, rec.seq.upper()))
    write_fasta(dirs.assemblies / f"{name}.fasta", out)
    return out


def _circular_count(records: list[Record]) -> int:
    return sum(1 for r in records if r.header.endswith(" circular=true"))


_CIRCULAR_TAGS = re.compile(r"(?:\bXO:i:1\b|circular[-=](?:yes|true|y)(?![a-z]))", re.I)


def count_circular(contigs: Path) -> int:
    """Contigs flagged circular in their header (myloasm: circular-yes; Raven: XO:i:1)."""
    return sum(1 for rec in read_records(contigs) if _CIRCULAR_TAGS.search(rec.header))


def _publish_graph(name: str, gfa: Path, dirs: AssemblyDirs, log_file: Path) -> None:
    if not gfa.is_file():
        return
    dirs.graphs.mkdir(parents=True, exist_ok=True)
    target = dirs.graphs / f"{name}.gfa"
    shutil.copyfile(gfa, target)
    if which("Bandage"):
        try:  # A picture of the graph is nice to have; its failure does not fail the sample
            run(["Bandage", "image", str(target), str(target.with_suffix(".png"))], log_file,
                env={"QT_QPA_PLATFORM": "offscreen"})
        except BaconError as exc:
            log.warning("Bandage could not draw the graph of %s: %s", name, str(exc).splitlines()[0])


def _assembly_stats(records: list[Record], reference_length: int) -> dict[str, object]:
    lengths = [len(r.seq) for r in records]
    total = sum(lengths)
    return {"Contigs": len(lengths), "Assembly_length": total, "Largest_contig": max(lengths),
            "Assembly_N50": n50(lengths), "Length_vs_reference": f"{total / reference_length:.3f}"}


def _flye_info(info: Path) -> dict[str, dict[str, str]]:
    rows = {}
    if info.is_file():
        with open(info) as fh:
            header = fh.readline().lstrip("#").strip().split("\t")
            for line in fh:
                values = dict(zip(header, line.rstrip("\n").split("\t")))
                if "seq_name" in values:
                    rows[values["seq_name"]] = values
    return rows


def assemble_flye(name: str, reads: Path, dirs: AssemblyDirs, log_dir: Path, *, genome_size: int,
                  read_type: str, min_overlap: int | None, iterations: int, threads: int,
                  reference_length: int) -> StepResult:
    work = dirs.work(name)
    shutil.rmtree(work, ignore_errors=True)
    dirs.root.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{name}.log"
    cmd = ["flye", f"--{read_type}", str(reads), "--genome-size", str(genome_size), "--threads", str(threads),
           "--out-dir", str(work), "--iterations", str(iterations)]
    if min_overlap:
        cmd += ["--min-overlap", str(min_overlap)]
    run(cmd, log_file, what=f"(assembling {name})")
    if not (work / "assembly.fasta").is_file():
        raise SampleFailed(f"Flye produced no assembly (see {work / 'flye.log'})")
    info = _flye_info(work / "assembly_info.txt")
    circular = {contig for contig, row in info.items() if row.get("circ.") == "Y"}
    records = _publish_assembly(name, work / "assembly.fasta", dirs, rename=False, circular=circular)
    _publish_graph(name, work / "assembly_graph.gfa", dirs, log_file)
    stats = _assembly_stats(records, reference_length)
    stats["Circular_contigs"] = _circular_count(records)
    covered = [(int(row["length"]), float(row["cov."])) for row in info.values()
               if row.get("length", "").isdigit() and _is_float(row.get("cov.", ""))]
    total = sum(length for length, _ in covered)
    stats["Assembly_depth"] = f"{sum(length * cov for length, cov in covered) / total:.1f}" if total else "NA"
    return StepResult(dirs.assemblies / f"{name}.fasta", stats)


def _is_float(text: str) -> bool:
    try:
        float(text)
    except ValueError:
        return False
    return True


def assemble_myloasm(name: str, reads: Path, dirs: AssemblyDirs, log_dir: Path, *, threads: int,
                     reference_length: int) -> StepResult:
    work = dirs.work(name)
    shutil.rmtree(work, ignore_errors=True)
    dirs.root.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{name}.log"
    run(["myloasm", str(reads), "--output-dir", str(work), "--threads", str(threads)], log_file,
        what=f"(assembling {name})")
    raw = work / "assembly_primary.fa"
    if not raw.is_file():
        raise SampleFailed(f"myloasm produced no assembly (see {log_file})")
    records = _publish_assembly(name, raw, dirs, rename=False)
    _publish_graph(name, work / "final_contig_graph.gfa", dirs, log_file)
    stats = _assembly_stats(records, reference_length)
    stats["Circular_contigs"] = _circular_count(records)
    return StepResult(dirs.assemblies / f"{name}.fasta", stats)


def assemble_samtools(name: str, reads: Path, reference: Path, dirs: AssemblyDirs, log_dir: Path, *,
                      threads: int, reference_length: int, fill_gaps: bool) -> StepResult:
    """Templated assembly: samtools consensus of the reads aligned to the reference with minimap2."""
    work = dirs.work(name)
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    log_file = log_dir / f"{name}.log"
    bam = work / "aligned.bam"
    run([["minimap2", "-ax", "map-ont", "--secondary=no", "-t", str(threads), str(reference), str(reads)],
         ["samtools", "sort", "-@", str(max(1, threads // 4)), "-o", str(bam), "-"]], log_file,
        what=f"(aligning {name})")
    run(["samtools", "index", str(bam)], log_file)
    consensus = work / "consensus.fasta"
    # -a: every reference position, N where no read; ONT R10.4 SUP error profile; insertions shown, deletions
    # removed, so the consensus follows the sample, not the reference.
    run(["samtools", "consensus", "-a", "-X", "r10.4_sup", "--show-ins", "yes", "--show-del", "no",
         "-d", "3", "-f", "fasta", "-o", str(consensus), str(bam)], log_file, what=f"(consensus {name})")
    records = [r for r in read_records(consensus) if r.seq]
    if fill_gaps:  # Replace N runs by the reference: only safe when the sample is known to match it there
        refs = {r.name: r.seq for r in read_records(reference)}
        records = [Record(r.header, _fill_from(r.seq, refs.get(r.name, ""))) for r in records]
        write_fasta(consensus, records)
    published = _publish_assembly(name, consensus, dirs, rename=False)
    stats = _assembly_stats(published, reference_length)
    n_bases = sum(r.seq.count("N") for r in published)
    stats["N_bases"] = n_bases
    notes = [f"{n_bases:,} N bases (no or ambiguous read support)"] if n_bases else []
    return StepResult(dirs.assemblies / f"{name}.fasta", stats, notes)


def _fill_from(seq: str, ref: str) -> str:
    """Leading/trailing N runs map to the reference ends; inner ones cannot be placed reliably, keep them."""
    if not ref:
        return seq
    core = seq.strip("N")
    if not core:  # No read support at all: nothing to anchor the reference on
        return seq
    lead = len(seq) - len(seq.lstrip("N"))
    trail = len(seq) - len(seq.rstrip("N"))
    return ref[:lead] + core + (ref[len(ref) - trail:] if trail else "")
