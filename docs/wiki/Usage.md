# Usage

```bash
bacon -r REFERENCE.fasta -i READS -o OUTPUT [options]
```

## Inputs

**Reference** (`-r`): a fasta file, gzipped or not, with one or more sequences, for example a complete
chloroplast genome. It is used to bait the reads, as the template of the templated assembly, and as the
anchor of the SNP comparison. A reference of the same species works best; a closely related species works for
baiting and de novo assembly. BACoN uses an upper-case copy (soft-masked, lower-case regions are compared like
the rest) in which IUPAC ambiguity codes are `N`; the same applies to the assemblies and to `--add-genomes`.

The reference can also be a GenBank file (`.gb`, `.gbk`, `.gbff`, `.genbank`, gzipped or not), such as NCBI's
record of a plastome: its sequence (`ORIGIN`) is the reference and its features annotate the report. Each
record is named after its `VERSION` (`NC_008096.2`), or its `LOCUS` name without one, followed by its
`DEFINITION`, as in NCBI's fasta of the same record: `reference.fasta` in the output folder is then identical to
the one made from that fasta, so switching between the two reruns nothing. `run_info.json` records the MD5 of
the file given, fasta or GenBank.

**Annotation** (`--annotation`): a GenBank or GFF3 file (`.gff`, `.gff3`; gzipped or not) describing the
reference, used in the report only: the genome map gets the genes, the regions of a plastome and the effect of
each SNP ([Outputs](Outputs#reporthtml)). It is not part of any checkpoint: adding, changing or removing it
reruns nothing. The default is the reference itself when it is a GenBank file. The annotated sequences must
have the names of the reference's sequences (a single annotated sequence of the same length as a single
reference sequence is accepted whatever its name); features on other sequences or beyond the end of a sequence
are ignored with a warning. BACoN copies the file, uncompressed, to `OUTPUT/annotation.gb` or
`OUTPUT/annotation.gff3`, so that the report can be rebuilt after the folder is moved.

**Reads**, as fastq or fasta, gzipped or not (`.fastq`, `.fq`, `.fasta`, `.fa`, `.fna`, `.fas`, with or
without `.gz`). Three ways to give them:

| Input | Samples |
|---|---|
| `-i sample.fastq.gz` | one sample, named after the file (`sample`) |
| `-i folder/` | each sequence file directly in the folder is a sample named after the file; each subfolder is a sample named after the subfolder, made of all the sequence files it contains, in any depth. MinKNOW's `fastq_pass/` works as is (`barcode01/`, `barcode02/`, ...); `unclassified/` and `mixed/` are skipped |
| `--sample-sheet samples.tsv` | a TSV or CSV file with the columns `sample` and `file`; several files per sample separated by `;` or on several rows; relative paths start from the sheet's folder; lines starting with `#` are ignored |

Sample names may contain letters, digits and `.` `_` `+` `-`. Two inputs giving the same sample name, or a
sample mixing fasta and fastq files, are errors. Symbolic links are a quick way to rename samples.

## Options

| Option | Default | Description |
|---|---|---|
| `-r`, `--reference` | | Reference fasta, or GenBank (its sequence is used and its features annotate the report); required |
| `-i`, `--input` | | Reads file or folder (this or `--sample-sheet`) |
| `--sample-sheet` | | TSV/CSV with `sample` and `file` columns |
| `-o`, `--output` | | Output folder (required) |
| `--annotation` | the reference, if GenBank | GenBank or GFF3 annotation of the reference, for the report only (genes, regions, SNP effects); never reruns a step |
| `-b`, `--baiting-method` | `minimap2` | `minimap2`: reads with an alignment to the reference; `bbduk`: reads sharing a k-mer (with `--hdist` mismatches) |
| `-k`, `--kmer-size` | 31 | BBDuk k-mer size (at most 31) |
| `--hdist` | 1 | BBDuk: mismatches allowed in a k-mer (0, 1 or 2). Each one multiplies BBDuk's memory: with 2, a 155 kb plastome needs about 14 GB per sample |
| `--keep-bam` | off | Keep the sorted BAM of the baited reads (minimap2) |
| `--min-read-length` | 500 | Shorter reads are discarded |
| `--keep-percent` | 95 | Filtlong keeps the best reads, up to this percentage of the bases (fastq only: fasta reads have no qualities and are selected by length) |
| `--target-depth` | 100 | Filtlong keeps at most this depth of the best reads (depth = bases / genome size, `-s` or the reference length) |
| `-a`, `--assembly-method` | `samtools` | `samtools` (templated), `flye` or `myloasm` (de novo); see [Methods](Methods) |
| `--template-gaps` | `n` | Templated assembly: reference positions covered by fewer than three reads are `N`; `reference` copies the reference into such positions at the sequence ends only |
| `--read-type` | `nano-hq` | Flye: `nano-hq` for R10/Q20 reads (<5% error; Flye's advice for R9 Guppy 5+ reads is also `nano-hq`), `nano-raw` for R9 reads basecalled with Guppy < 5, `nano-corr` for corrected reads |
| `--min-size` | automatic | Flye minimum read overlap |
| `-s`, `--size` | reference length | Expected genome size, for Flye and for `--target-depth` |
| `--snp-method`, `-snp` | `ska` | `ska` (SKA2 split k-mers), `parsnp` (core-genome alignment), `none` (stop after the assembly) |
| `--ska-min-freq` | 1.0 | SKA2: fraction of the genomes that must contain a variant's context; 1 = core SNPs, lower = pan-genome SNPs (like kSNP) |
| `--add-genomes` | | Finished genomes (fasta) to include in the comparison, such as published plastomes; named after their file |
| `--tree` | `fasttree` | `fasttree` (GTR, SH-like supports from 100 resamples) or `iqtree` (model selection, 1000 ultrafast bootstraps) |
| `--redo` | | Rerun this step and the following ones: `bait`, `filter`, `assemble`, `compare` |
| `-t`, `--threads` | all | Total threads, shared between the samples processed in parallel |
| `-p`, `--parallel` | 2 | Samples processed at the same time |
| `-m`, `--memory` | 85% of RAM | Total memory for BBDuk, in GB, divided between the samples processed in parallel (`-p`) |
| `--debug` | | Verbose log |

## Resuming and changing parameters

Each step records the parameters it ran with in `OUTPUT/.checkpoints/`. Running the same command again skips
every finished step; changing a parameter reruns that step and the ones after it, and only them. For example,
after a first run with the defaults:

```bash
bacon -r ref.fasta -i reads/ -o out/ -a flye              # reuses baiting and filtering, assembles with Flye
bacon -r ref.fasta -i reads/ -o out/ -a flye --snp-method parsnp   # reuses the Flye assemblies
```

Each comparison is written to its own folder (`4_compared/ska/`, `4_compared/parsnp/`,
`4_compared/ska_0.5/`), so several can be kept side by side. A sample that failed is retried on the next run;
samples added to the input are processed without redoing the others' baiting and filtering (their
assemblies are compared again). The files of a sample removed from the input stay in the output folder, unused;
after switching from a de novo to the templated assembly, the de novo assembly graphs of each sample are
removed.

`--redo STEP` forces a step to run again, for example after installing a newer assembler.

A run stopped with Ctrl-C (or SIGTERM, SIGHUP: `kill`, a closed terminal, a job scheduler's time limit) stops
the programs it started; each sample finished before the interruption is recorded, so the next run redoes only
the samples that were still running. A run started with `nohup` keeps running when the terminal is closed. Two
runs cannot use the same output folder at the same time: the second one stops with an error.

An output folder can be moved or copied and resumed from its new place. A copy that does not keep the files'
times (`cp -r` without `-a`, `rsync` without `-t`) keeps the baiting, but filtering, assembly and comparison
run again: BACoN recognizes each step's input by its size and time.

## Performance

The 28 potato samples of the [tutorial](Tutorial) (7 GB of whole-genome reads) take 3 minutes with the default
templated assembly (`-t 32 -p 8`): 1.5 minutes to bait, 1.2 minutes to filter, 17 seconds to assemble, and about
a second to compare. With Flye (`-t 48 -p 12`), the assembly takes 10 minutes: 2–4 minutes per sample, and 8
minutes for the slowest. Baiting reads every input read once; the other steps work on the baited reads only.
The largest single process used 0.7 GB with the templated assembly or Flye, and 2 GB with myloasm.
