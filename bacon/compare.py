"""Comparing assemblies: SNPs (SKA2 or Parsnp), SNP distances and a tree."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from bacon import BaconError
from bacon.newick import ladderize, midpoint_root, parse, to_newick, to_svg
from bacon.seqio import Record, read_records, write_fasta
from bacon.tools import run, which

log = logging.getLogger(__name__)

NUCLEOTIDES = frozenset("ACGT")


def clean_alignment(src: Path, dst: Path, rename: dict[str, str] | None = None) -> list[Record]:
    """Uppercase, and replace anything that is not A, C, G, T or '-' by N (like snippy-clean_full_aln)."""
    keep = frozenset("ACGT-")
    records = []
    for rec in read_records(src):
        seq = "".join(c if c in keep else "N" for c in rec.seq.upper())
        name = (rename or {}).get(rec.name, rec.name)
        records.append(Record(name, seq))
    if len({len(r.seq) for r in records}) > 1:
        raise BaconError(f"Sequences of {src} have different lengths: not an alignment")
    write_fasta(dst, records)
    return records


def snp_distances(records: list[Record]) -> tuple[list[str], list[list[int]]]:
    """Pairwise number of positions where both sequences have a nucleotide (A, C, G, T) and they differ."""
    records = sorted(records, key=lambda r: (r.name != "Reference", r.name))
    names = [r.name for r in records]
    n = len(records)
    matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a, b = records[i].seq, records[j].seq
            d = sum(1 for x, y in zip(a, b) if x != y and x in NUCLEOTIDES and y in NUCLEOTIDES)
            matrix[i][j] = matrix[j][i] = d
    return names, matrix


def write_distances(path: Path, names: list[str], matrix: list[list[int]]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(["snp-dists", *names]) + "\n")
        for name, row in zip(names, matrix):
            fh.write("\t".join([name, *map(str, row)]) + "\n")


# ---------------------------------------------------------------------------------------------------------------
# Parsnp
# ---------------------------------------------------------------------------------------------------------------

def _parsnp_names(assemblies: dict[str, Path], reference: Path) -> dict[str, str]:
    """Parsnp names sequences after their file (and adds '.ref' to the reference); map them back."""
    names = {path.name: name for name, path in assemblies.items()}
    names[reference.name + ".ref"] = "Reference"
    names[reference.name] = "Reference"
    return names


def run_parsnp(reference: Path, assemblies: dict[str, Path], out_dir: Path, log_dir: Path, *,
               threads: int) -> tuple[Path, Path, list[str]]:
    shutil.rmtree(out_dir, ignore_errors=True)
    log_file = log_dir / "parsnp.log"
    # -c: keep every genome (by default Parsnp silently drops genomes too distant from the reference).
    run(["parsnp", "-r", str(reference), "-d", *map(str, assemblies.values()), "-o", str(out_dir),
         "-p", str(threads), "-c"], log_file, what="(Parsnp)")
    xmfa, ggr = out_dir / "parsnp.xmfa", out_dir / "parsnp.ggr"
    if not xmfa.is_file() or not ggr.is_file():
        raise BaconError(f"Parsnp produced no alignment (assemblies too different from the reference?); "
                         f"see {log_file}")
    run(["harvesttools", "-x", str(xmfa), "-M", str(out_dir / "parsnp.core.raw.fasta")], log_file)
    run(["harvesttools", "-i", str(ggr), "-S", str(out_dir / "parsnp.snps.raw.fasta")], log_file)
    rename = _parsnp_names(assemblies, reference)
    clean_alignment(out_dir / "parsnp.core.raw.fasta", out_dir / "parsnp.core.fasta", rename)
    clean_alignment(out_dir / "parsnp.snps.raw.fasta", out_dir / "parsnp.snps.fasta", rename)
    for raw in ("parsnp.core.raw.fasta", "parsnp.snps.raw.fasta"):
        (out_dir / raw).unlink()
    return out_dir / "parsnp.core.fasta", out_dir / "parsnp.snps.fasta", []


# ---------------------------------------------------------------------------------------------------------------
# SKA2 (split k-mers, reference-free)
# ---------------------------------------------------------------------------------------------------------------

def wrap_circular(rec: Record, kmer: int, force: bool = False) -> Record:
    """Append the first k-1 bases of a circular sequence to its end (see run_ska)."""
    if (force or rec.header.endswith("circular=true")) and len(rec.seq) > kmer:
        return Record(rec.header, rec.seq + rec.seq[: kmer - 1])
    return rec


def run_ska(reference: Path, assemblies: dict[str, Path], out_dir: Path, log_dir: Path, *, threads: int,
            min_freq: float, kmer: int = 31) -> tuple[Path, Path, list[str]]:
    """Split k-mer alignment of the assemblies and the reference.

    `min_freq` is the fraction of genomes that must contain a split k-mer for its variant to be kept: 1 gives
    core SNPs, lower values keep SNPs missing from some genomes (a pan-genome alignment, gaps as '-').
    """
    shutil.rmtree(out_dir, ignore_errors=True)
    out_dir.mkdir(parents=True)
    log_file = log_dir / "ska.log"
    # Split k-mers spanning the start/end junction of a circular contig are lost, and with them any variant
    # within k bases of it: extend circular contigs by their first k-1 bases. The reference is treated as
    # circular when most assemblies are.
    inputs = out_dir / "inputs"
    inputs.mkdir()
    genomes: dict[str, Path] = {}
    circular_samples = 0
    for name, path in assemblies.items():
        records = list(read_records(path))
        if any(r.header.endswith("circular=true") for r in records):
            circular_samples += 1
        genomes[name] = inputs / f"{name}.fasta"
        write_fasta(genomes[name], [wrap_circular(r, kmer) for r in records])
    ref_records = list(read_records(reference))
    wrap_ref = circular_samples > len(assemblies) / 2
    genomes = {"Reference": inputs / "Reference.fasta", **genomes}
    write_fasta(genomes["Reference"], [wrap_circular(r, kmer, force=wrap_ref) for r in ref_records])
    table = out_dir / "input.tsv"
    table.write_text("".join(f"{name}\t{path}\n" for name, path in genomes.items()))
    run(["ska", "build", "-o", str(out_dir / "ska"), "-k", str(kmer), "-f", str(table), "--threads", str(threads)],
        log_file, what="(ska build)")
    raw = out_dir / "ska.raw.aln"
    run(["ska", "align", "--min-freq", f"{min_freq:g}", "--filter", "no-const", "-o", str(raw),
         "--threads", str(threads), str(out_dir / "ska.skf")], log_file, what="(ska align)")
    aln = out_dir / "ska.snps.fasta"
    clean_alignment(raw, aln)
    raw.unlink()
    shutil.rmtree(inputs)
    return aln, aln, []


# ---------------------------------------------------------------------------------------------------------------
# Trees
# ---------------------------------------------------------------------------------------------------------------

def build_tree(alignment: Path, out_dir: Path, log_dir: Path, *, method: str, threads: int) -> Path:
    """Build a tree from an alignment; write tree.nwk (midpoint-rooted) and tree.svg. Returns tree.nwk."""
    raw = out_dir / f"{method}.tree"
    if method == "fasttree":
        exe = which("FastTree") or "FastTree"
        run([exe, "-nt", "-gtr", "-boot", "100", str(alignment)], log_dir / "fasttree.log", stdout=raw,
            what="(FastTree)")
    elif method == "iqtree":
        exe = which("iqtree") or "iqtree"
        prefix = out_dir / "iqtree"
        run([exe, "-s", str(alignment), "-m", "MFP", "-B", "1000", "-T", str(threads), "--prefix", str(prefix),
             "-redo", "--seed", "12345"], log_dir / "iqtree.log", what="(IQ-TREE)")
        shutil.copyfile(prefix.with_suffix(".contree"), raw)
    else:
        raise ValueError(method)
    tree = midpoint_root(parse(raw.read_text()))
    ladderize(tree)
    out = out_dir / "tree.nwk"
    out.write_text(to_newick(tree) + "\n")
    label = "SH-like support (FastTree)" if method == "fasttree" else "ultrafast bootstrap (IQ-TREE)"
    (out_dir / "tree.svg").write_text(to_svg(tree, f"{alignment.name} — midpoint-rooted, {label}"))
    return out
