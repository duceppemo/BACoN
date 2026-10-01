# Changelog

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
- The VCF of a comparison made by BACoN 0.3.1 now uses the reference as SKA2 used it (extended when circular).
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
- A soft-masked or IUPAC reference rewrites `reference.fasta`, so every step runs again once.
- Comparisons with `--add-genomes` are redone once.
- `-b bbduk` baiting runs again once, for the new `--hdist` default.

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
