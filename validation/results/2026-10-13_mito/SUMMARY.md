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
| mitogenome, templated | **189–289** | **19,227–278,015** | 1, 367.5 kb | — | **29–34** | **5–8** |
| mitogenome, Flye | the plastome was assembled (154 kb; 1% of the truth) | — | 1–3 | **94–96%** | 28 | 0 |
| mitogenome, templated, `--filter-by-reference` | **0** | **105–162** | 1, 367.8 kb | — | **0** | **0** |
| mitogenome, Flye, `--filter-by-reference` | 0–6 | — | 3–4, 370–474 kb | 0–8.3% | 0–1 | 1 |

**Why the default fails on this mitogenome.** In one line (MA32-119), the baited reads covered the mitogenome
about 800x. Filtlong, which prefers long reads, kept 408 reads (median 83 kb) of which only 9% of the bases
aligned to the mitogenome: long reads with a stretch of mitochondrial-like sequence (nuclear copies; plastid reads,
on the plastid-derived stretches of the mitogenome). The estimated depth read 100x; the real depth was 3–14x.
With `--filter-by-reference` (Filtlong `--assembly`: reads ranked by their k-mers found in the reference), the
1,166 kept reads aligned in full, 100x.

**Homopolymers.** The templated plastomes have 15–33 one-base indels against the truth, all in runs of 8 or more
A or T, one base short; the reference itself differs from the truth by 1–10 indels. Substitutions: none. SKA2
compares substitutions only.

**Competitive baiting (tested, not kept).** Baiting the mitogenome against the plastome as a decoy (keeping the
reads whose alignments match more bases of the mitogenome) brought the templated mitogenomes to 128–20,130 `N`
and 0–21 substitutions, and Flye to the mitogenome plus 2–6 Mb of nuclear contigs. With `--filter-by-reference`,
adding the decoy changed nothing (the same `N` counts, the same substitutions); it was dropped.

## `--filter-by-reference` on the other data

Simulated scenarios (`run_validation.sh` with every bacon call given the flag): [`SCORES_*_filter.md`](.)

- plastid, linear, rdna: templated assemblies and every SNP distance unchanged. Flye assemblies changed by a few
  bases or contigs either way (s07: 0 to 5 errors; s09: 76 to 70 errors per 100 kb; l01 and l04, t01 and t04:
  completeness swapped), within Flye's run-to-run variation ([0.3.8 record](../2026-10-13_v0.3.8/SUMMARY.md)).
- numt: templated unchanged (`N` at the NUPT sites of the 25x samples: 721 and 29 by default, 749 and 27 with the
  flag); Flye better: no extra contig for the 10x sample, all 21 SNP distances exact (16 by default).
- heteroplasmy: Flye false SNPs at a 0.4 and 0.5 minor fraction 2 and 5 by default, 0 and 4 with the flag;
  templated unchanged.

Potato plastomes (the [Tutorial](../../../docs/wiki/Tutorial.md) data: 28 cultivars, R9.4.1, plastid-enriched
by baiting only): [`potato_filter.txt`](potato_filter.txt).

- Templated: 326 of 406 distances identical, the others within 2 SNPs; pairs of one lineage at a non-zero
  distance 69 of 166 (59 by default); the largest distance within a lineage 13 and the smallest between lineages
  66 in both. One cultivar (Argo) gained a block of 985 `N` at 52,420 (reads spanning sequence that differs from
  the reference rank low); median `N` 122 (133 by default).
- Flye: all 406 distances identical; 11 circular plastomes (17 by default).

Decision: an option, off by default, recommended for mitogenomes and whole-genome reads.

## `Off_reference_pct`

The share of a de novo assembly's bases that align nowhere on the reference (`minimap2 -x asm20`), with a note
above 10%. Computed on existing assemblies:

| Assemblies | n | Median | Max | Above 10% |
|---|---|---|---|---|
| simulated plastid, linear, rdna; Flye and myloasm (0.3.8 record) | 34 | 0.0% | 0.3% | 0 |
| simulated numt, Flye (two runs) | 12 | 0–10.4% | 49.6–55.1% | 5 (the 10x and 25x samples with extra contigs) |
| potato plastomes, Flye | 28 | 0.0% | 0.0% | 0 |
| potato plastomes, myloasm | 28 | 9.3% | 46.4% | 13 (extra contigs at 1–2x depth that align to the potato mitogenome) |
| potato mitogenomes (MN104801–3 as reference), Flye | 28 | 81.9% | 83.9% | 28 (the plastome, assembled from mitogenome baits) |
| Col-0 mitogenomes, Flye, default | 6 | 94.2% | 95.7% | 6 (the plastome) |
| Col-0 mitogenomes, Flye, `--filter-by-reference` | 6 | 0.0% | 8.3% | 0 |

## Heteroplasmy and nuclear copies (simulated)

[`SCORES_heteroplasmy.md`](SCORES_heteroplasmy.md), [`SCORES_numt.md`](SCORES_numt.md): scenarios described in
[`validation/README.md`](../../README.md), run with the defaults. Summarized in Methods (Heteroplasmy, Nuclear
copies):

- Heteroplasmy, 80x: the templated assembly gave the major allele at every site up to a minor fraction of 0.2,
  `N` from 0.3 (5 of 30 sites at 0.3, 22 at 0.5), never the minor allele; SNP distances exact. Flye gave the
  major allele up to 0.3; at 0.4 and 0.5, minor alleles or both alleles in different contigs, and 2 and 5 false
  SNPs.
- Nuclear copies (organelle 50x, nuclear DNA 0.5–25x): nearly every read covering a copy is baited. Templated:
  unaffected up to 10x; at 25x, `N` at the sites where the copies differ (721 of 1,703), never the nuclear base,
  and distances err low (2 SNPs for a truth of 5). Flye: extra contigs of the copies (assemblies 1.5–3.3 times
  the reference), some sites with both bases in different contigs.
