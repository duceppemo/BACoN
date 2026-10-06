# Changelog

## Unreleased

### Added
- Reference annotations: `--annotation FILE` (GenBank or GFF3, gzipped or not), or a GenBank file as `-r`
  (its sequence is the reference, named as in NCBI's fasta of the record, and its features the annotation).
  The report's genome map then shows the genes of each strand, the LSC/IRb/SSC/IRa regions of a plastome
  derived from its inverted repeats, and labels for the genes with the most SNPs, and a new sortable table
  gives each SNP's region, gene, context and effect on the coding sequence (codon and amino-acid change,
  synonymous/missense/nonsense/stop lost or retained/start lost or retained; translation table 11 unless the annotation
  says otherwise). The annotation is for the report only: it is not part of any checkpoint, so adding it to a
  finished run reruns nothing; it is copied to `OUTPUT/annotation.gb` or `.gff3` and recorded in
  `run_info.json`. `snps.vcf` is unchanged.
- Sample metadata: `--metadata FILE` (TSV or CSV with a `sample` column and any other columns), and the columns
  of a sample sheet other than `sample` and `file` (`--metadata` wins column by column). The report shows the
  columns in the samples table (sortable, numbers as numbers), and one column, `--color-by` or the first with at
  most 8 distinct values that is not free text (`none` for no colours), colours the tree (a circle and the
  value after each name) and the heatmap (a second band on both axes, with a legend), with a table of the groups
  of identical genomes against its values. The metadata columns are also in the MultiQC samples table. For the
  report only: not part of any checkpoint, copied to `OUTPUT/metadata.tsv` and recorded in `run_info.json`.

### Changed
- `report.html` has a new look (light and dark themes, numbered figures) and new figures: bar charts of the
  depth and of the `N` bases of each sample; the tree drawn with a coloured square for each group of identical
  genomes, the supports and a scale bar in SNPs; a heatmap of the distances with colour classes fitted to their
  range and the groups of identical genomes as coloured bands; and a genome map with the positions of the VCF's
  SNPs and, for templated assemblies, the `N` bases along the reference.

## 0.3.5 (2026-10-02)

### Changed
- BACoN stopped by SIGTERM or SIGHUP exits with 143 or 129 (128 + the signal), Ctrl-C with 130.
- SNP alignments (`ska.snps.fasta`, `parsnp.snps.fasta`) keep only the columns with at least two different
  nucleotides: with `--ska-min-freq` below 1, columns made of one nucleotide and gaps were counted as SNP sites
  (more than 95% of them in a simulated test) and given to the tree program; with any setting, a column whose
  only variation was an ambiguous base in one genome was counted too. The SNP distances do not change; the
  number of SNP sites and the trees' input do. Comparisons are not redone on resume: use `--redo compare`.
- The default `-t` and `-m` are the CPUs and memory BACoN may use under a job scheduler, `taskset` or in a
  container, not those of the whole machine.
- A templated assembly notes the reference sequences that no read covers (they are not in the assembly).

### Fixed
- VCF from Parsnp: positions where the reference's base is N (an ambiguity code in the reference) were written as
  a SNP of every genome.
- Midpoint rooting dropped the support of the split on which the root is placed (often the deepest one) from
  `tree.nwk`, `tree.svg` and the report.
- An output folder inside the input folder (`-i reads -o reads/bacon`) worked once, then every resume failed with
  an unrelated message (its files were taken for a sample); it is now refused.
- `-b bbduk`: a resumed run gave each sample its share of `-m` as if all samples were baited at the same time.
- A sample with no read baited lost its numbers of raw reads and bases in `summary.tsv`.
- Files in hidden folders inside a sample's folder were used.
- Read names with non-ASCII characters could not be matched, so their reads were not baited.
- A checkpoint of 0.3.3 or earlier with a sample named `distances` stopped the run.
- The report: the reference's file name was not escaped; the tree of a Parsnp comparison was said to be built on
  SNPs (it is built on the core-genome alignment); browsers that darken pages could make the tree invisible.
  MultiQC: the heatmap was said to be in tree order when there was no tree.
- `python -m bacon.report` on a folder that is not a BACoN output gave a traceback.
- `run()` left the `bacon` logger at the INFO level.
- `nohup bacon ... &` stopped when the terminal was closed (0.3.4 replaced nohup's ignored SIGHUP with its own
  handler); an ignored SIGHUP or SIGTERM now stays ignored.
- A resumed output folder of 0.3.3 or earlier could still reuse a result made from an earlier input (new reads
  baited, then a run interrupted before filtering): such results are now reused only if newer than their input,
  and get the input signatures of 0.3.4.
- The final error of a run (all samples failed, a failed comparison, "Interrupted") was not written to
  `bacon.log`.
- `bacon.cli.main` failed when called from a thread other than the main one, and left its signal handlers
  installed after returning.
- An interruption waited for the reads being copied or counted by BACoN itself (not by a program) to finish.
- A moved output folder of 0.3.3 or earlier whose path contained the name of a BACoN folder (such as
  `/data/1_extracted/out`) ran every step again.
- The upgrade notes: with `-b bbduk`, every step runs again once after an upgrade from 0.3.2, not only baiting.

## 0.3.4 (2026-10-01)

### Changed
- A VCF written by BACoN 0.3.2 or earlier is rewritten on resume, from the comparison's files (no need for
  `--redo compare`).
- SIGTERM and SIGHUP (`kill`, a closed terminal, a scheduler's time limit) stop BACoN like Ctrl-C: the programs
  it started are stopped. In 0.3.3 they kept running.

### Fixed
- 0.3.3 ran every step again on output folders of 0.3.2 (its checkpoints differed even with minimap2).
- A sample's result could be reused after its input changed, when a run was interrupted in between: new reads
  baited, then the run stopped before filtering; or a step run with other parameters, interrupted, and run
  again with the first ones. Each result now records the input it was made from, and so does the comparison.
- A moved or copied output folder resumed with the files of the original folder.
- After an interruption, the programs killed were reported as failures of their samples.
- On a file system without locks (some NFS or SMB mounts), BACoN said another run was using the folder.
- The report's methods said BBDuk allowed two mismatches, whatever `--hdist`.

### Upgrading to 0.3.4
- From 0.3.2: output folders resume, and their VCF is rewritten. Steps run again once only for a soft-masked or
  IUPAC reference (every step), `--add-genomes` (the comparison) and `-b bbduk` (every step, as baiting runs
  again for the new `--hdist` default).
- From 0.3.3: output folders baited with minimap2 run every step again once (0.3.3 changed the checkpoints by
  mistake); with `-b bbduk`, they resume.

## 0.3.3 (2026-10-01)

### Added
- `--hdist`: the mismatches BBDuk allows in a k-mer; the default is now 1 (was 2), which baits the same reads
  in the example with a fraction of the memory. BBDuk running out of memory gives a hint.

### Fixed
- An interrupted run (Ctrl-C) lost the samples finished in the step it was in; each sample is now recorded
  when it finishes. The programs still running, and their children (Flye's minimap2), are stopped.
- Two runs in the same output folder at the same time overwrote each other's files; the second one now stops.
- The traceback of an unexpected error in a sample is written to that sample's log.
- A file listed twice for one sample in a sample sheet doubled its reads; it is used once, with a warning.
- An `--add-genomes` genome could take the name of a sample that failed.
- SKA2: added genomes counted towards the rule that decides whether the genomes are circular, and were not
  extended when the reference was; the rule now counts the samples' assemblies, and the added genomes are
  extended with the reference.
- `--keep-bam` read all the read names through Python; `summary.tsv`, the distances, the VCF and the MultiQC
  files are written atomically; the report shows no inline heatmap above 150 genomes (the file is linked).
- A soft-masked (lower-case) reference lost the SNPs of its lower-case regions from `snps.vcf` (`ska map` and
  HarvestTools are case-sensitive), and IUPAC codes in a reference or added genome were read by SKA2 as fixed
  bases (false SNPs): the reference, the assemblies and the added genomes are now upper case with N for any
  ambiguity code.
- An unexpected error in one sample (a truncated or corrupt gzip file, a non-ASCII read name) stopped the whole
  run without a summary; it now fails that sample only, with a clear message. An unexpected error in the
  comparison still writes the summary, `run_info.json` and the report.
- `--add-genomes` made the comparison run again on every resume.
- Reads from fasta files were all kept when none reached `--min-read-length`.
- A failure in a pipe of programs was attributed to the upstream program killed by SIGPIPE.
- Parsnp: a warning when the core genome is less than half of the reference (a short genome can shrink it for
  every genome).
- Switching from a de novo to the templated assembly left the de novo assembly graphs.
- The VCF of a comparison made by BACoN 0.3.1 or 0.3.2 now uses the reference as SKA2 used it (extended when circular).
- Genotypes such as `0/1` or `1|0:35` are renumbered correctly when an `N` allele is removed.
- `example/run_example.sh` checks the SKA2 comparison, and its output folder is ignored by git.
- VCF: `N` (ambiguous or missing bases, such as the N of templated assemblies) was written as an alternate
  allele, and Parsnp wrote one record per `N`; it is now a missing genotype, and records left without an
  alternate allele are dropped, as are records where the reference's own base is ambiguous.
- VCF from SKA2: SNPs within 15 bp of the ends of circular genomes were in the alignment but not in the VCF;
  `ska map` now uses the reference as SKA2 used it (`ska_reference.fasta`), and positions are folded back.
- VCF from Parsnp: HarvestTools' undeclared INFO value `NA` made bcftools warn, and `bcftools norm` fail.
- A comparison made by BACoN 0.3.1, or whose VCF failed, never got a VCF on later runs; it now does, from the
  comparison's files, without redoing the comparison.
- A failed `ska map` left `snps.raw.vcf`; a malformed VCF from the tools stopped the run or gave an invalid
  file: it is now reported and the run goes on.

### Upgrading from 0.3.2
- Every step runs again once: the checkpoints changed (by mistake; 0.3.4 resumes folders of 0.3.2, with the
  exceptions below).
- A soft-masked or IUPAC reference rewrites `reference.fasta`, so every step runs again once.
- Comparisons with `--add-genomes` are redone once.
- `-b bbduk`: every step runs again once, as baiting runs again for the new `--hdist` default.

## 0.3.2 (2026-09-26)

### Added
- `snps.vcf` in the comparison folder: the SNPs of every genome relative to the reference (SKA2: `ska map`;
  Parsnp: HarvestTools), with contig lengths in the header, one column per genome
  ([#1](https://github.com/duceppemo/BACoN/issues/1)).

### Fixed
- MultiQC: the files of several runs in one search path were merged (samples pooled, only one heatmap kept);
  the sections are now named after the output folder.
- IQ-TREE renamed samples with `+` in their name, which then did not match in the report and MultiQC heatmaps;
  the tree programs now see placeholder names and the tree gets the sample names back.
- The report counted genomes as identical when they were connected through zero distances only (distances skip
  positions with `N`), and claimed identity when no SNP site was compared.
- The reference MD5 in `run_info.json` and the report was that of BACoN's copy, not of the file given.
- The methods paragraph of the report: fasta reads (no Filtlong, no quality cut), `--template-gaps
  reference`, Flye's genome size and overlap, added genomes, failed comparisons, IQ-TREE's consensus tree,
  Parsnp's core-genome tree, and the circular-contig extension only when there are circular contigs.
- `python -m bacon.report` on a moved output folder lost the tree and distances; incomplete `run_info.json`
  or an empty `summary.tsv` no longer crash it.
- Tabs or line breaks in an error message could shift the columns of `summary.tsv`.
- The report keeps its colours when printed, and sorts text columns naturally.

## 0.3.1 (2026-09-26)

### Added
- `report.html`: a self-contained HTML report written at the end of every run: overview, sortable sample table
  with flagged values, tree, SNP-distance heatmap in tree order, identical genomes, a methods paragraph with
  program versions, and provenance. `python -m bacon.report OUTPUT` rebuilds it.
- MultiQC custom content: `bacon_samples_mqc.json` (table), `bacon_reads_mqc.json` (bases kept, filtered out
  and off-target), `bacon_distances_mqc.json` (SNP-distance heatmap).
- `run_info.json` records the reference file, number of sequences, length and MD5.

## 0.3.0 (2026-09-25)

A rewrite. Every program was reconsidered: those no longer maintained were replaced after comparing the
candidates on simulated data of known truth and on real data ([Validation](https://github.com/duceppemo/BACoN/wiki/Validation)).

### Changed
- Installable package with a `bacon` command (`pip install .`); `python bacon.py` still works. BACoN itself
  uses the Python standard library only (pandas, psutil, pysam and ete3 are no longer needed).
- The default assembly is now **templated**: the consensus of the reads aligned to the reference, with
  `samtools consensus` (`-a samtools`). It replaces Rebaler and was the most accurate method in every test.
  Flye remains available (`-a flye`), and myloasm is added (`-a myloasm`).
- The default SNP method is now **SKA2** (`--snp-method ska`), reference-free, exact on all simulated data and
  counting a SNP in an inverted repeat once. `--ska-min-freq` below 1 gives pan-genome SNPs, the role of kSNP3
  in BACoN 0.1. Parsnp remains (`--snp-method parsnp`), now with `-c`, which turns off its filter that drops divergent genomes.
- Output folders: `1_extracted`, `2_filtered`, `3_assembled`, `4_compared/<method>/` (there is no trimming
  step any more). Each comparison method writes to its own folder.
- Resuming is automatic and parameter-aware: finished steps are skipped when their parameters and inputs did
  not change; changing a parameter reruns that step and the following ones; failed samples are retried;
  adding samples processes only them (and the comparison). `--redo STEP` forces a step. The `done_*` files
  are replaced by `.checkpoints/`.
- Trees are midpoint-rooted and drawn as SVG by BACoN (`tree.svg`), without Qt.
- Sample names are the file names without their extension (`_pass` is no longer removed); a folder of files
  (such as MinKNOW's `barcode01/`) is one sample.
- `--kmer-size` defaults to 31, the largest value BBDuk accepts (99 made BBDuk fail).

### Added
- `summary.tsv`: reads and bases at each step, baited share, depth, assembly statistics, circular contigs,
  status and notes (low depth, unexpected length, `N` bases, failures) for every sample.
- `snp_distances.tsv`: pairwise SNP distances.
- `run_info.json`: command, settings, versions of every program, samples and results; `bacon.log`; per-sample
  logs of every program in `logs/`.
- `--add-genomes`: finished genomes (such as published plastomes) included in the comparison.
- `--sample-sheet` (several files per sample), `--tree iqtree`, `--read-type`, `--template-gaps`,
  `--target-depth`, `--keep-percent`, `--min-read-length`, `--debug`.
- A failed sample no longer stops the run; every program's exit status is checked.
- Tests (pytest, stub programs), continuous integration, a simulated example with known SNPs
  (`example/`), a validation suite (`validation/`), a bioconda recipe, and a wiki maintained in `docs/wiki/`.

### Removed
- Porechop (unmaintained; reads are trimmed by the basecaller), Shasta, Rebaler, Snippy (cannot be installed
  with current assemblers), PhaME (its version check rejects samtools 1.10 to 1.29), RAxML and the
  `requirements.txt` environment. See [Methods](https://github.com/duceppemo/BACoN/wiki/Methods#choices-that-changed-in-03).

### Fixed
- With minimap2 baiting, the BAM files of other samples processed in parallel could be deleted before their
  reads were extracted.
- `--threads` and `--parallel` were not checked (the bounds test could never be true).
- Snippy ran with a changed working directory shared by parallel threads.
- PhaME added the reference to the assembly folder, which later runs took for a sample.
- Shasta needed `pigz`, which was not in the environment.
- Programs failed silently: their output was discarded and exit codes ignored.
- Reads given as fasta failed at filtering (Filtlong needs qualities); they are now filtered by length.

## 0.2 (2023-03-20)

- Snippy (default) and PhaME as SNP methods; kSNP3 removed.

## 0.1 (2022-08-30)

- First release: baiting (minimap2, BBDuk), Porechop, Filtlong, Flye/Shasta/Rebaler, Parsnp and kSNP3.
