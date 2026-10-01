"""Comparing assemblies: SNPs (SKA2 or Parsnp), SNP distances and a tree."""

from __future__ import annotations

import logging
import re
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
               threads: int) -> tuple[Path, Path]:
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
    reference_length = sum(len(r.seq) for r in read_records(reference))
    core_length = next((len(r.seq) for r in read_records(out_dir / "parsnp.core.fasta")), 0)
    if core_length < 0.5 * reference_length:
        log.warning("Parsnp's core genome is %s bp, %.0f%% of the reference: an incomplete assembly limits the "
                    "comparison of every genome to that core; check the assemblies, or use --snp-method ska",
                    f"{core_length:,}", 100 * core_length / reference_length)
    return out_dir / "parsnp.core.fasta", out_dir / "parsnp.snps.fasta"


# ---------------------------------------------------------------------------------------------------------------
# SKA2 (split k-mers, reference-free)
# ---------------------------------------------------------------------------------------------------------------

def wrap_circular(rec: Record, kmer: int, force: bool = False) -> Record:
    """Append the first k-1 bases of a circular sequence to its end (see run_ska)."""
    if (force or rec.header.endswith("circular=true")) and len(rec.seq) > kmer:
        return Record(rec.header, rec.seq + rec.seq[: kmer - 1])
    return rec


SKA_REFERENCE = "ska_reference.fasta"


def write_ska_reference(reference: Path, assemblies: dict[str, Path], out: Path, kmer: int = 31) -> bool:
    """The reference as SKA2 uses it: extended by its first k-1 bases when most assemblies are circular.
    Returns whether it was extended."""
    circular = sum(1 for path in assemblies.values()
                   if any(r.header.endswith("circular=true") for r in read_records(path)))
    wrap = circular > len(assemblies) / 2
    write_fasta(out, [wrap_circular(r, kmer, force=wrap) for r in read_records(reference)])
    return wrap


def run_ska(reference: Path, assemblies: dict[str, Path], out_dir: Path, log_dir: Path, *, threads: int,
            min_freq: float, kmer: int = 31) -> tuple[Path, Path]:
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
    for name, path in assemblies.items():
        genomes[name] = inputs / f"{name}.fasta"
        write_fasta(genomes[name], [wrap_circular(r, kmer) for r in read_records(path)])
    genomes = {"Reference": inputs / "Reference.fasta", **genomes}
    write_ska_reference(reference, assemblies, genomes["Reference"], kmer)
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
    # The reference as SKA2 saw it (extended when circular): write_vcf maps to it, so that the VCF has the SNPs
    # near the ends that the alignment has.
    genomes["Reference"].replace(out_dir / SKA_REFERENCE)
    shutil.rmtree(inputs)
    return aln, aln


# ---------------------------------------------------------------------------------------------------------------
# VCF
# ---------------------------------------------------------------------------------------------------------------

def _genotype_fix(value: str, remap: dict[str, str]) -> str:
    """Renumber the alleles of a genotype field (GT first; '0', '1/2', '0|1' and extra ':' fields kept)."""
    gt, sep, rest = value.partition(":")
    alleles = re.split(r"([/|])", gt)
    return "".join(remap.get(a, a) if a not in "/|" else a for a in alleles) + sep + rest


def clean_vcf(raw: Path, out: Path, *, rename: dict[str, str], reference: Path, source: str) -> int:
    """Rewrite a VCF from SKA2 or HarvestTools into a plain SNP VCF. Returns the number of records.

    - Sample columns are named after the genomes and sorted; the reference's own column is removed.
    - The header has the contig names and lengths of `reference`, the source, and the GT format.
    - `N` is not an allele: genotypes pointing to it become missing (`.`) and it is removed from ALT; a record
      left without an alternate allele is dropped (positions a genome lacks, such as a deletion, and ambiguous
      bases). Records where the reference itself is not `0` (its own split k-mer is ambiguous) are dropped.
    - Positions beyond the end of a contig (the reference was extended by the start of a circular sequence)
      are folded back onto the start; duplicates are dropped.
    - HarvestTools' undeclared INFO value `NA` becomes `.`; its `N` filter is removed with the `N` allele.
    """
    lengths = {rec.name: len(rec.seq) for rec in read_records(reference)}
    meta: list[str] = []
    header: list[str] | None = None
    keep: list[int] = []
    ref_col: int | None = None
    records: dict[tuple[str, int], list[str]] = {}
    with open(raw) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith("##"):
                if not line.startswith(("##contig", "##source", "##fileformat")):
                    meta.append(line)
                continue
            fields = line.split("\t")
            if line.startswith("#CHROM"):
                names = [rename.get(n, n) for n in fields[9:]]
                ref_cols = [i for i, n in enumerate(names) if n == "Reference"]
                if len(ref_cols) != 1:
                    raise BaconError(f"{raw}: expected one reference column, found {len(ref_cols)} ({names})")
                ref_col = ref_cols[0]
                keep = sorted((i for i, n in enumerate(names) if n != "Reference"), key=lambda i: names[i])
                header = fields[:9] + [names[i] for i in keep]
                continue
            if header is None or ref_col is None:
                raise BaconError(f"{raw}: variant records before the #CHROM header line")
            if len(fields) < 9 + len(names):
                raise BaconError(f"{raw}: truncated record: {line[:80]}")
            genotypes = fields[9:]
            if genotypes[ref_col] != "0":
                continue
            alts = fields[4].split(",")
            remap, kept_alts = {}, []
            for index, allele in enumerate(alts, 1):
                if allele in ("N", ".", ""):
                    remap[str(index)] = "."
                else:
                    kept_alts.append(allele)
                    remap[str(index)] = str(len(kept_alts))
            if not kept_alts:
                continue
            chrom, pos = fields[0], int(fields[1])
            length = lengths.get(chrom)
            if length and pos > length:
                pos -= length
            if (chrom, pos) in records:
                continue
            filters = [f for f in fields[6].split(";") if f not in ("N", "")]
            row = [chrom, str(pos), fields[2], fields[3], ",".join(kept_alts), fields[5],
                   ";".join(filters) if filters else ("PASS" if fields[6] not in (".", "") else "."),
                   "." if fields[7] in ("NA", "") else fields[7], fields[8]]
            records[(chrom, pos)] = row + [_genotype_fix(genotypes[i], remap) for i in keep]
    if header is None:
        raise BaconError(f"{raw}: no #CHROM header line")
    if not any(m.startswith("##FORMAT=<ID=GT,") for m in meta):
        meta.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    order = {name: i for i, name in enumerate(lengths)}
    tmp = out.with_name(out.name + ".tmp")
    with open(tmp, "w") as dst:
        dst.write("##fileformat=VCFv4.2\n")
        dst.write(f"##source={source}\n")
        dst.writelines(f"##contig=<ID={name},length={length}>\n" for name, length in lengths.items())
        dst.writelines(m + "\n" for m in meta)
        dst.write("\t".join(header) + "\n")
        for key in sorted(records, key=lambda k: (order.get(k[0], len(order)), k[0], k[1])):
            dst.write("\t".join(records[key]) + "\n")
    tmp.replace(out)
    return len(records)


def write_vcf(method: str, reference: Path, out_dir: Path, log_dir: Path, *, threads: int,
              assemblies: dict[str, Path], source: str) -> tuple[Path, int]:
    """SNPs of every genome relative to the reference, as VCF (`snps.vcf` in the comparison folder)."""
    raw = out_dir / "snps.raw.vcf"
    out = out_dir / "snps.vcf"
    out.unlink(missing_ok=True)
    try:
        if method == "ska":
            mapped_to = out_dir / SKA_REFERENCE
            if not mapped_to.exists():  # A comparison made by BACoN 0.3.1: rebuild it by the same rule
                write_ska_reference(reference, assemblies, mapped_to)
            run(["ska", "map", str(mapped_to), str(out_dir / "ska.skf"),
                 "-f", "vcf", "-o", str(raw), "--threads", str(threads)], log_dir / "ska.log", what="(ska map)")
            rename: dict[str, str] = {}
        else:
            run(["harvesttools", "-i", str(out_dir / "parsnp.ggr"), "-V", str(raw)], log_dir / "parsnp.log",
                what="(HarvestTools VCF)")
            rename = _parsnp_names(assemblies, reference)
        count = clean_vcf(raw, out, rename=rename, reference=reference, source=source)
    finally:
        raw.unlink(missing_ok=True)
    return out, count


# ---------------------------------------------------------------------------------------------------------------
# Trees
# ---------------------------------------------------------------------------------------------------------------

def build_tree(alignment: Path, out_dir: Path, log_dir: Path, *, method: str, threads: int) -> Path:
    """Build a tree from an alignment; write tree.nwk (midpoint-rooted) and tree.svg. Returns tree.nwk.

    The tree programs see placeholder names (IQ-TREE rewrites characters such as '+'), which are replaced by the
    genome names in tree.nwk and tree.svg.
    """
    records = list(read_records(alignment))
    names = {f"g{i:05d}": rec.name for i, rec in enumerate(records, 1)}
    safe = out_dir / "tree_input.fasta"
    write_fasta(safe, [Record(key, rec.seq) for key, rec in zip(names, records)])
    raw = out_dir / f"{method}.tree"
    try:
        if method == "fasttree":
            exe = which("FastTree") or "FastTree"
            run([exe, "-nt", "-gtr", "-boot", "100", str(safe)], log_dir / "fasttree.log", stdout=raw,
                what="(FastTree)")
        elif method == "iqtree":
            exe = which("iqtree") or "iqtree"
            prefix = out_dir / "iqtree"
            run([exe, "-s", str(safe), "-m", "MFP", "-B", "1000", "-T", str(threads), "--prefix", str(prefix),
                 "-redo", "--seed", "12345"], log_dir / "iqtree.log", what="(IQ-TREE)")
            shutil.copyfile(prefix.with_suffix(".contree"), raw)
        else:
            raise ValueError(method)
    finally:
        safe.unlink(missing_ok=True)
    tree = parse(raw.read_text())
    for leaf in tree.leaves():
        if leaf.name not in names:
            raise BaconError(f"{method} returned an unexpected sequence name {leaf.name!r}")
        leaf.name = names[leaf.name]
    raw.write_text(to_newick(tree) + "\n")  # The program's tree, with the genome names
    tree = midpoint_root(tree)
    ladderize(tree)
    out = out_dir / "tree.nwk"
    out.write_text(to_newick(tree) + "\n")
    label = "SH-like support (FastTree)" if method == "fasttree" else "ultrafast bootstrap (IQ-TREE)"
    (out_dir / "tree.svg").write_text(to_svg(tree, f"{alignment.name} — midpoint-rooted, {label}"))
    return out
