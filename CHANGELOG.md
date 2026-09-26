# Changelog

## Unreleased

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
