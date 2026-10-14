# Mitogenomes, nuclear copies, heteroplasmy; `--filter-by-reference` and `Off_reference_pct` (2026-10-13)

The validation of 0.3.0 and 0.3.8 tested plastomes, linear targets and tandem arrays only. This record adds:
real mitogenomes (*Arabidopsis* Col-0, with a truth), two simulated scenarios (heteroplasmy, nuclear copies of
organelle DNA), and the two changes made from what they showed: `--filter-by-reference` (an option) and the
off-reference share of de novo assemblies (`Off_reference_pct`).

BACoN: the development version after 0.3.9 (the changes of this record). Programs: those of 0.3.8
([`../2026-10-13_v0.3.8/versions.txt`](../2026-10-13_v0.3.8/versions.txt): minimap2 2.31, Filtlong 0.3.1, Flye
2.9.6, myloasm 0.7.0, samtools 1.24, SKA2 0.5.1). Truth: TIPPo 2.4 (with Flye 2.9.5), MUMmer 4.0.1. One machine:
2 Intel Xeon Gold 6226R (64 threads), 503 GB of RAM, Linux.

## *Arabidopsis* Col-0, plastome and mitogenome

Data: six mutation-accumulation lines of Col-0 (MA32-59, -69, -79, -89, -99, -119; ENA PRJEB112076), each
sequenced with ONT R10.4.1 (whole-genome reads, 11–19 Gb per line) and PacBio HiFi. Lines of one ecotype a few
generations apart: their organelle genomes should be identical, or nearly.

Truth: each line's organelles assembled from its HiFi reads by TIPPo (reference-free; for MA32-69, whose
mitochondrial assembly failed in TIPPo, Flye on TIPPo's mitochondrial reads). BACoN ran with the references
NC_000932.1 (plastome) and NC_037304.1 (mitogenome, 367,808 bp), with the defaults and with
`--filter-by-reference`, templated and Flye. Each assembly was aligned to its line's truth with nucmer
(`delta-filter -1`, `show-snps`): substitutions between ACGT bases, indel events, and the share of the truth
covered. The mitogenome truths are 4 contigs of 413–446 kb (Flye keeps some repeat copies apart), so a complete
mitogenome covers about 0.8–0.94 of them. Scripts: [`scripts/`](scripts/); all values:
[`col0_scores.tsv`](col0_scores.tsv).

| Run | Substitutions vs truth | `N` bases | Contigs, length | Off-reference | SNPs to the reference | Max SNPs between lines |
|---|---|---|---|---|---|---|
| plastome, templated | 0 | 0 | 1, 154.4 kb | — | 0 | 0 |
| plastome, Flye | 0 | — | 1–4 | 0% | 0 | 0 |
| plastome, templated, `--filter-by-reference` | 0 | 0 | 1, 154.4 kb | — | 0 | 0 |
| plastome, Flye, `--filter-by-reference` | 0 | — | 2, 154.4 kb (truth covered in full by all six) | 0% | 0 | 0 |
| mitogenome, templated | **189–289** | **19,227–278,015** | 1, 367.5–368.4 kb | — | **29–34** | **5–8** |
| mitogenome, Flye | the plastome was assembled (154 kb; 1% of the truth) | — | 1–3 | **94–96%** | 28 | 0 |
| mitogenome, templated, `--filter-by-reference` | **0** | **105–162** | 1, 367.8 kb | — | **0** | **0** |
| mitogenome, Flye, `--filter-by-reference` | 0–6 | — | 3–4, 370–474 kb | 0–8.3% | 0–1 | 1 |

**Why the default fails on this mitogenome** ([`col0_mt_read_filter.tsv`](col0_mt_read_filter.tsv)). In one line
(MA32-119), the baited reads covered the mitogenome 819x (median). Filtlong, which prefers long reads, kept 391
reads (median 90 kb) of which 5.7% of the bases aligned to the mitogenome: long reads with a stretch of
mitochondrial-like sequence (nuclear copies; plastid reads, on the plastid-derived stretches of the mitogenome).
The estimated depth read 100x; the median depth of the kept reads on the mitogenome was 3x. With
`--filter-by-reference` (Filtlong `--assembly`: reads ranked by their k-mers found in the reference), the 1,166
kept reads were 99.7% aligned, 100x.

**Homopolymers** ([`col0_plastome_indels.tsv`](col0_plastome_indels.tsv)). The six default templated plastomes
have 116 indels against the truth (15–23 each; 15–33 in the other plastome runs, `col0_scores.tsv`): 114 are
deletions of one base (96) or two (18) in runs of 11–16 A or T, the runs being one or two bases short; the two
others are insertions in runs of 5. The reference itself differs from the truths by 1–10 one-base indels (16 in
all). Substitutions: none. SKA2 compares substitutions only.

**Competitive baiting (tested, not kept).** Baiting the mitogenome against the plastome as a decoy (keeping the
reads whose alignments match more bases of the mitogenome) brought the templated mitogenomes to 128–20,130 `N`
and 0–21 substitutions, and Flye to the mitogenome plus 2–6 Mb of nuclear contigs. With `--filter-by-reference`,
adding the decoy changed nothing (the same `N` counts, the same substitutions); it was dropped.

## `--filter-by-reference` on the other data

Simulated scenarios (`run_validation.sh` with every bacon call given the flag): [`SCORES_*_filter.md`](.)

- plastid, linear, rdna: templated assemblies and every SNP distance unchanged. Flye assemblies changed by a few
  bases or contigs either way (s07: 0 to 5 errors; s09: 76 to 70 errors per 100 kb; l01 0.9999 to 0.9819 and l04
  0.9593 to 0.9999 complete; t01 0.27 to 0.14 and t04 0.27 to 0.41), within Flye's run-to-run variation ([0.3.8 record](../2026-10-13_v0.3.8/SUMMARY.md)).
- numt: templated unchanged (`N` at the NUPT sites of the 25x samples: 721 and 29 by default, 749 and 27 with the
  flag); Flye better: all 21 SNP distances exact (16 by default); the 10x sample of genome A
  went from 3 contigs (185 kb) to 2 (149 kb), and the 25x samples from 397 and 377 kb to 339 and 184 kb.
- heteroplasmy: Flye false SNPs at a 0.4 and 0.5 minor fraction 2 and 5 by default, 0 and 4 with the flag;
  templated unchanged.

Potato plastomes (the [Tutorial](../../../docs/wiki/Tutorial.md) data: 28 cultivars, R9.4.1, plastid-enriched
by baiting only): [`potato_filter.txt`](potato_filter.txt).

- Templated: 326 of 406 distances identical, the others within 2 SNPs; pairs of one lineage at a non-zero
  distance 69 of 166 (59 by default); the largest distance within a lineage 13 and the smallest between lineages
  66 in both. One cultivar (Argo) gained a run of 985 `N` at 52,420 (reads spanning sequence that differs from
  the reference rank low; by default, 151 `N` there in short runs); median `N` 122 (133 by default).
- Flye: all 406 distances identical; 11 circular plastomes (17 by default).

Decision: an option, off by default, recommended for mitogenomes and whole-genome reads.

## `Off_reference_pct`

The share of a de novo assembly's bases that align nowhere on the reference (`minimap2 -x asm20`), with a note
above 10%. Computed on existing assemblies ([`off_reference.tsv`](off_reference.tsv), every value):

| Assemblies | n | Median | Max | Above 10% |
|---|---|---|---|---|
| simulated plastid, linear, rdna; Flye and myloasm (0.3.8 record) | 34 | 0.0% | 0.3% | 0 |
| simulated numt, Flye | 6 | 10.4% | 49.6% | 3 (n03, n04, n06: the samples with extra contigs) |
| simulated numt, Flye, `--filter-by-reference` | 6 | 4.8% | 47.6% | 2 (n04, n06); n03, with an extra 29 kb contig, at 9.7% |
| potato plastomes, Flye | 28 | 0.0% | 0.0% | 0 |
| potato plastomes, myloasm | 28 | 9.3% | 46.4% | 13 (extra contigs at 1–2x depth that align to the potato mitogenome) |
| potato mitogenomes (MN104801–3 as reference), Flye | 28 | 81.9% | 83.9% | 28 (the plastome, assembled from mitogenome baits) |
| Col-0 mitogenomes, Flye, default | 6 | 94.2% | 95.7% | 6 (the plastome) |
| Col-0 mitogenomes, Flye, `--filter-by-reference` | 6 | 0.0% | 8.3% | 0 |
| Col-0 plastomes, Flye, with and without `--filter-by-reference` | 12 | 0.0% | 0.0% | 0 |

Assemblies of the target alone stayed at 0–0.3% where there was nothing else to assemble (the simulated
plastid, linear and rdna scenarios, the potato and Col-0 plastomes with Flye). The Col-0 mitogenomes with
`--filter-by-reference` reached 2.7% and 8.3%; the share of a mitogenome that a de novo assembler resolves
differently from the reference is larger. The note at 10% missed one assembly with an extra contig (n03 above).

## Heteroplasmy and nuclear copies (simulated)

[`SCORES_heteroplasmy.md`](SCORES_heteroplasmy.md), [`SCORES_numt.md`](SCORES_numt.md): scenarios described in
[`validation/README.md`](../../README.md), run with the defaults. Summarized in Methods (Heteroplasmy, Nuclear
copies):

- Heteroplasmy, 80x: the templated assembly gave the major allele at every site up to a minor fraction of 0.2,
  `N` from 0.3 (5 of 30 sites at 0.3, 22 at 0.5), never the minor allele; SNP distances exact. Flye gave the
  major allele up to 0.3; at 0.4 and 0.5, minor alleles or both alleles in different contigs, and 2 and 5 false
  SNPs.
- Nuclear copies (organelle 50x, nuclear DNA 0.5–25x): nearly every read covering a copy is baited. Templated:
  unaffected up to 2x; at 10x, one `N`, which cost one SNP of a distance (4 for a truth of 5); at 25x, `N` at the sites where the copies differ (721 of 1,703), never the nuclear base,
  and distances err low (2 SNPs for a truth of 5). Flye: extra contigs of the copies (assemblies 1.5 times
  the reference at 10x, 3.1–3.3 at 25x), some sites with both bases in different contigs.
