# Usage

```bash
bacon -r REFERENCE.fasta -i READS -o OUTPUT [options]
```

## Inputs

**Reference** (`-r`): a fasta file, gzipped or not, with one or more sequences, for example a complete
chloroplast genome. It is used to bait the reads, as the template of the templated assembly, and as the
anchor of the SNP comparison. A reference of the same species works best; a closely related species works for
baiting and de novo assembly.

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
| `-r`, `--reference` | | Reference fasta (required) |
| `-i`, `--input` | | Reads file or folder (this or `--sample-sheet`) |
| `--sample-sheet` | | TSV/CSV with `sample` and `file` columns |
| `-o`, `--output` | | Output folder (required) |
| `-b`, `--baiting-method` | `minimap2` | `minimap2`: reads with an alignment to the reference; `bbduk`: reads sharing a k-mer (up to 2 mismatches) |
| `-k`, `--kmer-size` | 31 | BBDuk k-mer size (at most 31) |
| `--keep-bam` | off | Keep the sorted BAM of the baited reads (minimap2) |
| `--min-read-length` | 500 | Shorter reads are discarded |
| `--keep-percent` | 95 | Filtlong keeps this percentage of the best reads (fastq only: fasta reads have no qualities and are selected by length) |
| `--target-depth` | 100 | Filtlong keeps at most this depth of the best reads (depth = bases / genome size, `-s` or the reference length) |
| `-a`, `--assembly-method` | `samtools` | `samtools` (templated), `flye` or `myloasm` (de novo); see [Methods](Methods) |
| `--template-gaps` | `n` | Templated assembly: reference positions covered by fewer than three reads are `N`; `reference` copies the reference into such positions at the sequence ends only |
| `--read-type` | `nano-hq` | Flye: `nano-hq` for R10/Q20 reads (<3% error; Flye's advice for R9 Guppy 5+ reads is also `nano-hq`), `nano-raw` for R9 reads basecalled with Guppy < 5, `nano-corr` for corrected reads |
| `--min-size` | automatic | Flye minimum read overlap |
| `-s`, `--size` | reference length | Expected genome size, for Flye and for `--target-depth` |
| `--snp-method`, `-snp` | `ska` | `ska` (SKA2 split k-mers), `parsnp` (core-genome alignment), `none` (stop after the assembly) |
| `--ska-min-freq` | 1.0 | SKA2: fraction of the genomes that must contain a variant's context; 1 = core SNPs, lower = pan-genome SNPs (like kSNP) |
| `--add-genomes` | | Finished genomes (fasta) to include in the comparison, such as published plastomes; named after their file |
| `--tree` | `fasttree` | `fasttree` (GTR, SH-like supports from 100 resamples) or `iqtree` (model selection, 1000 ultrafast bootstraps) |
| `--redo` | | Rerun this step and the following ones: `bait`, `filter`, `assemble`, `compare` |
| `-t`, `--threads` | all | Total threads, shared between the samples processed in parallel |
| `-p`, `--parallel` | 2 | Samples processed at the same time |
| `-m`, `--memory` | 85% of RAM | Memory for BBDuk, in GB |
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
assemblies are compared again).

`--redo STEP` forces a step to run again, for example after installing a newer assembler.

## Performance

The 28 potato samples of the [tutorial](Tutorial) (7 GB of whole-genome reads) take 3 minutes with the default
templated assembly (`-t 32 -p 8`): 1.5 minutes to bait, 1.2 minutes to filter, 17 seconds to assemble, and about
a second to compare. With Flye (`-t 48 -p 12`), the assembly takes 10 minutes: 2–4 minutes per sample, and 8
minutes for the slowest. Baiting reads every input read once; the other steps work on the baited reads only.
The largest single process used 0.7 GB with the templated assembly or Flye, and 2 GB with myloasm.
