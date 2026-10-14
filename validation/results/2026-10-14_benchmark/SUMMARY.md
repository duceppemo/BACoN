# BACoN and other organelle pipelines on *Arabidopsis* Col-0, ONT R10.4.1 (2026-10-14)

The same six Col-0 lines and HiFi truths as the [mitogenome record](../2026-10-13_mito/SUMMARY.md), assembled by
four other pipelines that take Nanopore reads: ptGAUL, TIPPo, Oatk and CLAW. Each pipeline's final assemblies were
scored against each line's truth as in that record, and the six lines compared with SKA2 as BACoN compares
genomes. BACoN: the version of commit 7c98957 (after 0.3.9, with `--filter-by-reference`); its assemblies are those
of the mitogenome record, rescored.

| Pipeline | Version | Organelles | Reference | Assembler | Options (besides threads) |
|---|---|---|---|---|---|
| BACoN | 7c98957 | one per run | baits the reads, and is the template (templated) | samtools consensus or Flye 2.9.6 | defaults; mitogenome with `--filter-by-reference` |
| [ptGAUL](https://github.com/Bean061/ptgaul) | 1.0.5 (bioconda) | plastome | baits the reads | Flye 2.7 (pinned by its recipe) | `-g 154478` |
| [TIPPo](https://github.com/Wenfei-Xian/TIPPo) | 2.4 | plastome, mitogenome | none | Flye 2.9.5 | `-p onthq -g organelle` |
| [Oatk](https://github.com/c-zhou/oatk) | 1.0 (bioconda) | plastome, mitogenome | none (HMMs of OatkDB v20230921, embryophyta) | its own (syncasm) | `-k 501 -c 30` |
| [CLAW](https://github.com/aaronphillips7493/CLAW) | git clone 2026-10-13, Snakemake 7.32.4 | plastome | baits the reads | Flye 2.9.6 and Unicycler 0.5.1 on 3,000 reads of ≥ 5 kb | `--nano-hq`, `chloroplast_size: 154478` |

Oatk is designed for HiFi reads: with its default k-mer size (1001), the trial on MA32-99 gave the plastome and no
mitogenome; k 501 gave both, and was used for all six. Not run: pipelines for short reads (GetOrganelle,
NOVOPlasty) or HiFi reads only (MitoHiFi).

Scoring ([`scripts/score.py`](scripts/score.py)): nucmer `--maxmatch` of each assembly against its line's truth;
substitutions (between ACGT bases) and indel events from the one-to-one alignments (`delta-filter -1`), the share of
the truth covered by alignments of at least 1 kb (`delta-filter -l 1000`; see the erratum below). A plastome
covered in full scores 1.0 whether its inverted repeat is in one copy or two. Distances
([`scripts/distances.py`](scripts/distances.py)): SKA2 core SNPs (k 31) of each pipeline's six assemblies and the
reference, BACoN's default comparison. All values: [`scores.tsv`](scores.tsv), [`distances.tsv`](distances.tsv)
(matrices in [`distances/`](distances/)), [`times.tsv`](times.tsv).

One machine: 2 Intel Xeon Gold 6226R (64 threads), 503 GB of RAM, Linux.

## Plastome

| Pipeline | Complete (truth covered ≥ 0.99) | Contigs, length | Substitutions vs truth | Indels vs truth | SNPs to the reference | Max SNPs between lines |
|---|---|---|---|---|---|---|
| BACoN, templated | 6 | 1, 154.4 kb | 0 | 15–23 | 0 | 0 |
| BACoN, Flye | 6 | 1–4 (three with one inverted-repeat copy, 128 kb) | 0 | 9–28 | 0 | 0 |
| ptGAUL | 6 | 1, 154.4 kb | 0 | 18–30 | 0 | 0 |
| TIPPo | 6 | 2 (the two orientations of the small single-copy region), 2 × 154.4 kb | 0 | 17–51 | 0 | 0 |
| Oatk | 4 | 1, 154.4 kb; MA32-59 139 contigs (0.63 covered), MA32-69 184 contigs (0.95) | 0–1 | 28–45 | 0 | 0 |
| CLAW, Flye | 5 | 1, 154.4 kb; MA32-99 one contig of 43.5 kb (0.29) | 0 (2 in the 43.5 kb contig) | 17–26 (5) | 0 | 0 |
| CLAW, Unicycler | 5 | 1, 154.5 kb (two lines); 2–5, 482–682 kb (three); MA32-99 one 11.6 kb contig that aligns nowhere on the truth | 0–5; MA32-69 176 | 15–33; MA32-69 222 | 0 | 0 |

The complete plastomes of BACoN, ptGAUL, TIPPo, Oatk and CLAW-Flye had no substitution against their truth;
CLAW-Unicycler's had 0–5, and 176 in MA32-69. Every pipeline gave the six lines identical plastomes (no SNP).
In BACoN's templated plastomes the indels are nearly all homopolymer runs one or two bases short (the
[mitogenome record](../2026-10-13_mito/SUMMARY.md)); the other pipelines' indels were not classified. SKA2
compares substitutions only. This version of CLAW keeps every contig of each assembler's
output (its size and identity checks are commented out of its Snakefile), so its Unicycler assemblies above 154 kb
hold the plastome and other contigs; its coverage tracks (deepTools `bamCoverage`) failed for every line on a
matplotlib incompatibility in its environment, after the assemblies were written.

## Mitogenome

| Pipeline | Truth covered | Contigs, length | Substitutions vs truth | Indels vs truth | SNPs to the reference | Max SNPs between lines |
|---|---|---|---|---|---|---|
| BACoN, templated, `--filter-by-reference` | 0.999–1.0 | 1, 367.8 kb | 0 | 4–8 | 0 | 0 |
| BACoN, Flye, `--filter-by-reference` | 0.999–1.0 | 3–4, 370–474 kb | 0–6 | 4–25 | 0–1 | 1 |
| BACoN, templated (default) | 0.19–0.93 (`N`) | 1, 367.5–368.4 kb | 189–289 | 50–64 | 29–34 | 8 |
| BACoN, Flye (default) | 0 (the plastome) | 1–3, 128–154 kb | — | — | — | — |
| TIPPo | 0.97–1.0, 4 of 6 lines (Flye failed on MA32-59 and -69: "No disjointigs were assembled") | 4–5, 380–475 kb | 0–18 | 5–29 | 0–2 | 1 |
| Oatk | 0.11–0.70 | 1–10, 48–238 kb | 0–7 | 2–12 | 0 | 1 |

The mitogenome truths are 4 contigs of 413–446 kb (HiFi Flye keeps some repeat copies apart); BACoN's templated
mitogenome follows the reference's single 367.8 kb sequence. Oatk's mitochondrial contigs at k 501 are partial on
every line; its distances cover only the part common to the six assemblies.

## Wall time and memory

Each pipeline alone on the machine with 16 threads ([`times.tsv`](times.tsv)). BACoN: the six lines one after
the other (`-p 1`), fresh runs ([`scripts/run_bacon_timing.sh`](scripts/run_bacon_timing.sh)), the mean per line.
The others: one line, MA32-119 (12.0 GB of gzipped reads; the six lines' mean 12.7 GB), one pipeline after the
other ([`scripts/run_solo.sh`](scripts/run_solo.sh)).

| Pipeline | Organelles | Wall time per line | Peak memory |
|---|---|---|---|
| BACoN, templated | plastome | 9.0 min | 1.3 GB |
| BACoN, Flye | plastome | 10.3 min | 1.3 GB |
| BACoN, templated, `--filter-by-reference` | mitogenome | 10.4 min | 1.3 GB |
| BACoN, Flye, `--filter-by-reference` | mitogenome | 14.7 min | 1.3 GB |
| ptGAUL | plastome | 14.1 min | 7.6 GB |
| CLAW (Flye and Unicycler) | plastome | 30.4 min | 4.4 GB |
| Oatk, k 501 | both | 45.9 min | 39.8 GB |
| TIPPo | both | 58.8 min | 84.6 GB |

BACoN builds one organelle per run: both, templated, take 19.4 min per line. Most of BACoN's time is the baiting
of the whole-genome reads: in the templated plastome run, 42 of the 54 min for the six lines (minimap2), then 12
min of Filtlong and under a minute for the consensus and SKA2. CLAW's time includes its coverage-track steps, which
failed (above), and its `snakemake` exited with an error after writing the assemblies, as in the six-line run.

The six-line runs of ptGAUL (32 threads), TIPPo (three lines at a time), CLAW and the Oatk trial overlapped on
the machine, and Oatk ran three lines at a time: their times in [`times.tsv`](times.tsv) are for a shared machine
(ptGAUL 15–55 min per line, TIPPo 60–130, Oatk 35–107, CLAW 155 min for the six lines), and peak memory per line
reached 13.4 GB (ptGAUL), 119.5 GB (TIPPo) and 54.3 GB (Oatk).

**BACoN's fresh runs against the record's.** The templated assemblies were identical (plastome and mitogenome).
Flye differs from run to run: the fresh plastomes had the same scores (0 substitutions, the same indels) and other
contig breaks in four lines (MA32-99: 5 contigs instead of 4); the fresh mitogenomes, 0–6 substitutions again
(0, 1, 0, 0, 6, 0 instead of 5, 3, 0, 0, 6, 0), 1–5 contigs of 367–441 kb, 0.999–1.0 of the truth covered. The
tables above give the record's runs.

## Erratum: `truth_covered` in the 2026-10-13_mito record

The mitogenome record computed the share of the truth covered from the one-to-one alignments
(`delta-filter -1`). These keep a single copy of each repeat, so a complete genome scored low: four complete
plastomes at 0.83 (one inverted-repeat copy dropped), one at 0.595, and the complete mitogenomes with
`--filter-by-reference` at 0.80–0.94 (templated) and 0.82–0.99 (Flye). The record explained the latter by the
truths' extra repeat copies ("a complete mitogenome covers about 0.8–0.94 of them"); that explanation was wrong.
Rescored with every alignment of at least 1 kb: every complete plastome and every `--filter-by-reference`
mitogenome covers 0.998–1.0 of its truth. The default templated mitogenomes move from 0.21–0.87 to 0.19–0.93
and the default Flye "mitogenomes" (the plastome) from 0.008–0.009 to 0 (the plastid-derived stretches of the mitogenome
are shorter than 1 kb or split). Substitutions, indels and every other value of that record are unchanged (same
alignments). Both values for every sample: [`erratum_truth_covered.tsv`](erratum_truth_covered.tsv).

## Files

- [`scripts/`](scripts/): `run_bacon_timing.sh`, `run_ptgaul.sh`, `run_tippo_ont.sh`, `run_oatk_trial.sh`,
  `run_oatk.sh`, `run_claw.sh`, `run_solo.sh` (the runs; `BENCH` is the benchmark folder: reads in
  `col0_lines/ont/`, references in `refs/`, truths in `col0/truth_tippo/`, made by the mitogenome record's
  scripts); `collect.sh` (the final assemblies of each pipeline), `score_all.sh` (`score.py`, `distances.py`),
  `summarize.py` (the tables).
