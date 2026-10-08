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
reruns nothing. The LSC/IRb/SSC/IRa band of a plastome comes from the annotated inverted repeats or, when the
annotation does not mark them (older RefSeq plastomes), from the large inverted repeat BACoN detects in the
reference sequence itself, when the regions have the layout of a plastome (the repeats at least 5% of the
sequence, the larger single-copy region at most 200 kb, the record not a mitochondrion or a chromosome); the
band is drawn with or without an annotation, the gene track only with one
([Methods](Methods#5-snp-effects)). The default is the reference itself when it is a GenBank file. The annotated sequences must
have the names of the reference's sequences (a single annotated sequence of the same length as a single
reference sequence is accepted whatever its name); features on other sequences or beyond the end of a sequence
are ignored with a warning, and an annotated sequence whose length (`LOCUS`, `##sequence-region`) is not the
reference's gets a warning asking whether it is the annotation of this reference. BACoN copies the file,
uncompressed, to `OUTPUT/annotation.gb` or `OUTPUT/annotation.gff3`, so that the report can be rebuilt after the
folder is moved. What is read from the annotation, and how SNP effects are derived, is in
[Methods](Methods#5-snp-effects).

**Sample metadata** (`--metadata`): a TSV or CSV file with a `sample` column (any case) naming the samples and
any other columns, such as a group, a cultivar, a site or a year, used in the report only
([Outputs](Outputs#reporthtml)): the columns are added to the samples table, and one of them colours the bar
charts, the tree and the heatmap (each sample has the same colour and shape everywhere in the report), with a
table of the groups of identical genomes against its values. Like the annotation, it is
not part of any checkpoint: adding, changing or removing it reruns nothing, and the report is rebuilt. The
columns of a sample sheet other than `sample` and `file` are metadata too (a sample on several rows must have the
same values; otherwise the first is kept, with a warning); when both are given, a column of `--metadata` replaces
the sheet's column of the same name (whatever its case), and its columns come first. Names and values are
stripped, and runs of whitespace inside a metadata value become one space (a sample sheet's file paths are kept
as written); empty cells, `NA`, `na` and `-` are missing values. Two columns of the same name, whatever their
case, are an error; a column named like one of the samples table's own columns (`Status`, `Note`, `Depth`...)
is shown as `Status (metadata)`. A sample listed twice keeps its first row (warning); rows naming no sample of
the run are reported; samples without a row get blank cells. Genomes given with `--add-genomes` may have a row
too, under their file name. BACoN writes the merged table, limited to the run's samples and added genomes, to
`OUTPUT/metadata.tsv`, so that the report can be rebuilt after the folder is moved; a later run without metadata
removes that copy (a `metadata.tsv` that BACoN did not write is kept, and the report uses it, saying that it was
found in the folder and not given to the run; `run_info.json` and the MultiQC table do not have it).

The file format is the same for a metadata file and a sample sheet, with one difference: blank lines are
ignored, and lines starting with `#` before the header are comments; after the header, a line starting with `#`
is a comment in a sample sheet too (a sample is left out by commenting its line) but a row in a metadata file
(so a value such as `#FF0000` may come first). A tab in the header makes the file a TSV: a tab always separates
cells, a cell entirely in quotes loses them (`"abc"` is `abc`, and `""` inside such a cell is one quote, as a
spreadsheet exports them) and other quotes are ordinary characters (`5" tube`); otherwise it is a CSV, whose
quoted values may contain commas but not line breaks (an unclosed quote is an error naming the line, not a value
that silently swallows the following rows).

```
sample      group      year  comment
alpha       A          2021  Identical to the reference by design
beta        A          2021  Carries six planted SNPs
gamma       B          2022  Shares beta's six SNPs, plus four of its own
delta       NA         2023  Ten SNPs of its own; group unknown
```

`--color-by COLUMN` chooses the column that colours every genome in the report (tree, heatmap, bar charts,
tables); `--color-by none` leaves them uncoloured (and needs no metadata). By
default it is the first column that can be coloured: a column with at most 48 distinct values that does not
look like free text, that is, whose distinct values are at most 30 characters long on average and, once 10 or
more samples have a value, do not outnumber half of them (so 28 samples can have up to 14 values, and 48 values
need 96 samples). Numbers with few distinct values (years, batches) are categories like any other. A requested
column that does not exist stops BACoN before any step; one that cannot be coloured gives a warning and
uncoloured figures. The values are marked in sorted order (numerically when they are all numbers), each with a
colour and a shape of its own: 12 colours from the colourblind-safe palettes of Okabe & Ito and of Paul Tol,
chosen and ordered to stay as distinct as possible under deuteranopia, protanopia and tritanopia (the first ones,
which a column with few values gets, the most distinct), and seven shapes (circle, triangle, square, diamond,
inverted triangle, plus, cross). The first 12 values take the colours in turn with the first four shapes in turn,
so that the values sharing a shape are four colours apart and stay distinct to a colourblind reader; the next 36
were chosen by a search to keep the values sharing a shape as distinct as possible. Up to 25 values, values
whose colours look alike to a colourblind reader therefore differ in shape; with more, a few pairs sharing a
shape may look alike (more so towards 48 values), and the legend and the value written after each name in the
tree and in the tables tell them apart. Consecutive values differ in both colour and shape. The markers follow
the sorted values, so adding or removing a value (a new sample, an edited sheet) can change the colour and shape
of the values after it. A genome without a value gets a hollow circle.

**Reads**, as fastq or fasta, gzipped or not (`.fastq`, `.fq`, `.fasta`, `.fa`, `.fna`, `.fas`, with or
without `.gz`). Three ways to give them:

| Input | Samples |
|---|---|
| `-i sample.fastq.gz` | one sample, named after the file (`sample`) |
| `-i folder/` | each sequence file directly in the folder is a sample named after the file; each subfolder is a sample named after the subfolder, made of all the sequence files it contains, in any depth (linked subfolders included; a subfolder without any is skipped, with a message). MinKNOW's `fastq_pass/` works as is (`barcode01/`, `barcode02/`, ...); `unclassified/` and `mixed/` are skipped |
| `--sample-sheet samples.tsv` | a TSV or CSV file with the columns `sample` and `file`; several files per sample separated by `;` or on several rows; relative paths start from the sheet's folder; blank lines and lines starting with `#` are ignored; any other column is sample metadata (above) |

Sample names may contain letters, digits and `.` `_` `+` `-`. Two inputs giving the same sample name, or a
sample mixing fasta and fastq files, are errors. Symbolic links are a quick way to rename samples.

## Options

| Option | Default | Description |
|---|---|---|
| `-r`, `--reference` | | Reference fasta, or GenBank (its sequence is used and its features annotate the report); required |
| `-i`, `--input` | | Reads file or folder (this or `--sample-sheet`) |
| `--sample-sheet` | | TSV/CSV with `sample` and `file` columns; other columns are metadata |
| `-o`, `--output` | | Output folder (required) |
| `--annotation` | the reference, if GenBank | GenBank or GFF3 annotation of the reference, for the report only (genes, regions, SNP effects); never reruns a step |
| `--metadata` | | TSV/CSV with a `sample` column and any other columns, for the report only (samples table, colours of every genome in the report); never reruns a step |
| `--color-by` | first usable column | Colours every genome in the report (tree, heatmap, bar charts, tables) by this metadata column (at most 48 distinct values, each with a colour and a shape; not free text), or `none` |
| `-b`, `--baiting-method` | `minimap2` | `minimap2`: reads with an alignment to the reference; `bbduk`: reads sharing a k-mer (with `--hdist` mismatches) |
| `-k`, `--kmer-size` | 31 | BBDuk k-mer size (at most 31) |
| `--hdist` | 1 | BBDuk: mismatches allowed in a k-mer (0, 1 or 2). Each one multiplies BBDuk's memory: with 2, a 155 kb plastome needs about 14 GB per sample (ignored, with a warning, with `-b minimap2`) |
| `--keep-bam` | off | Keep the sorted BAM of the baited reads (minimap2) |
| `--min-read-length` | 500 | Shorter reads are discarded |
| `--keep-percent` | 95 | Filtlong keeps the best reads, up to this percentage of the bases (fastq only: fasta reads have no qualities and are selected by length) |
| `--target-depth` | 100 | Filtlong keeps at most this depth of the best reads (depth = bases / genome size, `-s` or the reference length) |
| `-a`, `--assembly-method` | `samtools` | `samtools` (templated), `flye` or `myloasm` (de novo); see [Methods](Methods) |
| `--template-gaps` | `n` | Templated assembly: what to put where fewer than three reads cover the reference. `n`: an `N` (no base is guessed); `reference`: the reference's bases, at the ends of each sequence only (gaps inside stay `N`) |
| `--read-type` | `nano-hq` | Flye: `nano-hq` for R10/Q20 reads (<5% error; Flye's advice for R9 Guppy 5+ reads is also `nano-hq`), `nano-raw` for R9 reads basecalled with Guppy < 5, `nano-corr` for corrected reads |
| `--min-size` | automatic | Flye minimum read overlap |
| `-s`, `--size` | reference length | Expected genome size, for Flye and for `--target-depth` |
| `--snp-method`, `-snp` | `ska` | `ska` (SKA2 split k-mers), `parsnp` (core-genome alignment), `none` (stop after the assembly) |
| `--ska-min-freq` | 1.0 | SKA2: fraction of the genomes that must contain a variant's context; 1 = core SNPs, lower = pan-genome SNPs (like kSNP); see [Methods](Methods#4-comparison) |
| `--add-genomes` | | Finished genomes (fasta) to include in the comparison, such as published plastomes; named after their file |
| `--tree` | `fasttree` | `fasttree` (GTR, SH-like supports from 100 resamples) or `iqtree` (model selection, 1000 ultrafast bootstraps) |
| `--redo` | | Rerun this step and the following ones: `bait`, `filter`, `assemble`, `compare` |
| `-t`, `--threads` | all | Total threads, shared between the samples processed in parallel |
| `-p`, `--parallel` | 2 | Samples processed at the same time |
| `-m`, `--memory` | 85% of RAM | Total memory for BBDuk, in GB, divided between the samples processed in parallel (`-p`) |
| `--debug` | | Verbose log |
| `-v`, `--version` | | Print BACoN's version and exit |
| `-h`, `--help` | | Print the options and exit |

## Resuming and changing parameters

Each step records the parameters it ran with in `OUTPUT/.checkpoints/`. Running the same command again skips
every finished step; changing a parameter reruns that step and the ones after it, and only them (Flye's options
with another assembler rerun nothing, and `--keep-bam` only adds BAMs: it reruns only the samples without one).
For example, after a first run with the defaults:

```bash
bacon -r ref.fasta -i reads/ -o out/ -a flye              # reuses baiting and filtering, assembles with Flye
bacon -r ref.fasta -i reads/ -o out/ -a flye --snp-method parsnp   # reuses the Flye assemblies
```

Each comparison is written to its own folder (`4_compared/ska/`, `4_compared/parsnp/`,
`4_compared/ska_0.5/`), so several can be kept side by side. A sample that failed is retried on the next run;
samples added to the input are processed without redoing the others' baiting and filtering (their
assemblies are compared again). The files of a sample removed from the input (or whose line of the sample
sheet is commented out) stay in the output folder, unused, and are reused if it comes back unchanged; after
switching from a de novo to the templated assembly, the de novo assembly graphs of each sample are
removed.

`--redo STEP` forces a step to run again, for example after installing a newer assembler.

A run stopped with Ctrl-C (or SIGTERM, SIGHUP: `kill`, a closed terminal, a job scheduler's time limit) stops
the programs it started; each sample finished before the interruption is recorded, so the next run redoes only
the samples that were still running. A run started with `nohup` keeps running when the terminal is closed. Two
runs cannot use the same output folder at the same time: the second one stops with an error.

A run killed with SIGKILL (`kill -9`, the kernel's out-of-memory killer) cannot stop the programs it started:
each runs in its own process group, and they keep running (and writing to the output folder) until they end.
Their command lines name the output folder: list them with `pgrep -af /full/path/to/out` (or
`ps -ef | grep /full/path/to/out`) and, once the list holds only them, stop them with `pkill -f /full/path/to/out`
before running BACoN again on that folder.

An input file in the output folder that BACoN would overwrite (`-r`, `--annotation`, `--metadata`,
`--sample-sheet` named as one of the files BACoN writes there, e.g. `reference.fasta` or `metadata.tsv`) is an
error, unless it is already exactly the copy BACoN would write (`-r out/reference.fasta` of an earlier run).

An output folder can be moved or copied and resumed from its new place, whether or not the copy kept the files'
times (`cp -r`, `scp`, `rsync` without `-t`, an archive): BACoN recognizes each step's input by its size and
time in the folder where it was made, and by its size alone in a moved or copied folder.

## Performance

The 28 potato samples of the [tutorial](Tutorial) (7 GB of whole-genome reads) take 3 minutes with the default
templated assembly (`-t 32 -p 8`): 1.5 minutes to bait, 1.2 minutes to filter, 17 seconds to assemble, and about
a second to compare. With Flye (`-t 48 -p 12`), the assembly takes 10 minutes: 2–4 minutes per sample, and 8
minutes for the slowest. Baiting reads every input read once; the other steps work on the baited reads only.
The largest single process used 0.7 GB with the templated assembly or Flye, and 2 GB with myloasm.
