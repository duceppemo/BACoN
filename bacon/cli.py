"""Command-line interface for BACoN."""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

from bacon import BaconError, __version__
from bacon.pipeline import STEPS, Settings, default_memory_gb, run

log = logging.getLogger("bacon")


def _positive_int(value: str) -> int:
    try:
        i = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not an integer: {value!r}") from None
    if i < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1, got {value}")
    return i


def _percent(value: str) -> float:
    try:
        f = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a number: {value!r}") from None
    if not 0 < f <= 100:
        raise argparse.ArgumentTypeError(f"must be above 0 and at most 100, got {value}")
    return f


def _fraction(value: str) -> float:
    try:
        f = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a number: {value!r}") from None
    if not 0 < f <= 1:
        raise argparse.ArgumentTypeError(f"must be above 0 and at most 1, got {value}")
    return f


def _kmer(value: str) -> int:
    k = _positive_int(value)
    if k > 31:
        raise argparse.ArgumentTypeError(f"BBDuk k-mers are at most 31, got {value}")
    return k


def build_parser() -> argparse.ArgumentParser:
    max_cpu = os.cpu_count() or 1
    max_mem = default_memory_gb()
    parser = argparse.ArgumentParser(
        prog="bacon",
        description="Bait, Assemble and Compare Nanopore reads matching a reference sequence "
                    "(for example an organelle genome).",
    )
    io = parser.add_argument_group("input/output")
    io.add_argument("-r", "--reference", metavar="FILE.fasta", type=Path, required=True,
                    help="Reference sequence(s) used to bait the reads and to call SNPs (fasta, gzipped or not).")
    io.add_argument("-i", "--input", metavar="PATH", type=Path,
                    help="A fastq/fasta file, or a folder: each file is a sample named after the file, and each "
                         "subfolder is a sample named after the subfolder with all the files it contains (e.g. "
                         "MinKNOW's fastq_pass/barcode01/). Gzipped or not.")
    io.add_argument("--sample-sheet", metavar="FILE.tsv", type=Path,
                    help="Instead of -i: a TSV/CSV with the columns 'sample' and 'file' (several files per sample "
                         "separated by ';' or on several rows). Relative paths start from the sheet's folder.")
    io.add_argument("-o", "--output", metavar="DIR", type=Path, required=True,
                    help="Output folder. Rerunning into the same folder resumes: finished steps whose "
                         "parameters did not change are skipped.")

    bait = parser.add_argument_group("baiting")
    bait.add_argument("-b", "--baiting-method", choices=["minimap2", "bbduk"], default="minimap2",
                      help="minimap2 keeps reads that align to the reference; bbduk keeps reads sharing a k-mer "
                           "with it. Default: %(default)s")
    bait.add_argument("-k", "--kmer-size", metavar="31", type=_kmer, default=31,
                      help="K-mer size for bbduk (at most 31). Default: %(default)s")
    bait.add_argument("--keep-bam", action="store_true",
                      help="Keep the sorted BAM of the reads aligned to the reference (minimap2 only).")

    reads = parser.add_argument_group("read preparation")
    reads.add_argument("--min-read-length", metavar="500", type=_positive_int, default=500,
                       help="Discard reads shorter than this. Default: %(default)s")
    reads.add_argument("--keep-percent", metavar="95", type=_percent, default=95.0,
                       help="Keep this percentage of the best reads (Filtlong). Default: %(default)s")
    reads.add_argument("--target-depth", metavar="100", type=_positive_int, default=100,
                       help="Keep at most this depth of the best reads (Filtlong). Default: %(default)s")

    asm = parser.add_argument_group("assembly")
    asm.add_argument("-a", "--assembly-method", choices=["samtools", "flye", "myloasm"], default="samtools",
                     help="samtools: templated assembly, the consensus of the reads aligned to the reference "
                          "(most accurate for SNPs, fastest); flye or myloasm: de novo assembly (independent of "
                          "the reference, shows structure; myloasm keeps linear tandem arrays linear). "
                          "Default: %(default)s")
    asm.add_argument("--template-gaps", choices=["n", "reference"], default="n",
                     help="samtools: reference positions that no read covers are N, or at the sequence ends only, "
                          "copied from the reference. Default: %(default)s")
    asm.add_argument("--read-type", choices=["nano-hq", "nano-raw", "nano-corr"], default="nano-hq",
                     help="Flye read type: nano-hq for Guppy5+/Dorado reads (<5%% error), nano-raw for older "
                          "reads. Default: %(default)s")
    asm.add_argument("--min-size", metavar="BP", type=_positive_int, default=None,
                     help="Flye minimum read overlap. Default: automatic.")
    asm.add_argument("-s", "--size", metavar="BP", type=_positive_int, default=None,
                     help="Expected genome size, for Flye and read filtering. Default: the reference length.")

    cmp = parser.add_argument_group("comparison")
    cmp.add_argument("--snp-method", "-snp", choices=["ska", "parsnp", "none"], default="ska",
                     help="How to find SNPs between the assemblies: ska (SKA2 split k-mers, reference-free; core "
                          "or pan-genome SNPs, see --ska-min-freq) or parsnp (core-genome alignment to the "
                          "reference); 'none' stops after the assembly. Default: %(default)s")
    cmp.add_argument("--ska-min-freq", metavar="1.0", type=_fraction, default=1.0,
                     help="SKA2: fraction of genomes that must share a variant's context for it to be kept. "
                          "1 = core SNPs; lower values give a pan-genome alignment. Default: %(default)s")
    cmp.add_argument("--add-genomes", metavar="FASTA", type=Path, nargs="+", default=[],
                     help="Finished genomes (fasta) to include in the comparison, such as published plastomes; "
                          "each is named after its file.")
    cmp.add_argument("--tree", choices=["fasttree", "iqtree"], default="fasttree",
                     help="Tree builder: FastTree (GTR, SH-like support) or IQ-TREE (model selection, 1000 "
                          "ultrafast bootstraps). Default: %(default)s")

    run_group = parser.add_argument_group("run control")
    run_group.add_argument("--redo", choices=STEPS, default=None,
                           help="Run this step and the following ones again, even if already done.")
    run_group.add_argument("-t", "--threads", metavar=str(max_cpu), type=_positive_int, default=max_cpu,
                           help="Total number of threads, shared between parallel samples. Default: all "
                                "(%(default)s)")
    run_group.add_argument("-p", "--parallel", metavar="2", type=_positive_int, default=2,
                           help="Number of samples processed at the same time. Default: %(default)s")
    run_group.add_argument("-m", "--memory", metavar="GB", type=_positive_int, default=max_mem,
                           help="Memory in GB, for bbduk. Default: 85%% of the total (%(default)s)")
    parser.add_argument("--debug", action="store_true", help="Verbose logging.")
    parser.add_argument("-v", "--version", action="version", version=f"BACoN {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if (args.input is None) == (args.sample_sheet is None):
        parser.error("give exactly one of -i/--input or --sample-sheet")
    if args.keep_bam and args.baiting_method != "minimap2":
        parser.error("--keep-bam only applies to --baiting-method minimap2")
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    max_cpu = os.cpu_count() or 1
    threads = args.threads
    if threads > max_cpu:
        log.warning("Requested %d threads but only %d CPU(s) are available; using %d", threads, max_cpu, max_cpu)
        threads = max_cpu
    if args.parallel > threads:
        log.warning("--parallel (%d) exceeds --threads (%d); each sample gets 1 thread", args.parallel, threads)

    settings = Settings(
        reference=args.reference,
        output=args.output,
        input=args.input,
        sample_sheet=args.sample_sheet,
        baiting=args.baiting_method,
        kmer=args.kmer_size,
        keep_bam=args.keep_bam,
        min_read_length=args.min_read_length,
        keep_percent=args.keep_percent,
        target_depth=args.target_depth,
        assembler=args.assembly_method,
        read_type=args.read_type,
        min_size=args.min_size,
        template_gaps=args.template_gaps,
        genome_size=args.size,
        snp_method=args.snp_method,
        ska_min_freq=args.ska_min_freq,
        add_genomes=args.add_genomes,
        tree=args.tree,
        threads=threads,
        parallel=args.parallel,
        memory_gb=args.memory,
        redo=args.redo,
        command_line=list(sys.argv if argv is None else ["bacon", *argv]),
    )
    try:
        return run(settings)
    except BaconError as exc:
        log.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        log.error("Interrupted")
        return 130


if __name__ == "__main__":
    sys.exit(main())
