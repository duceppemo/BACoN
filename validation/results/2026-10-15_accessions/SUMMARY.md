# BACoN and other pipelines on 24 *Arabidopsis* accessions, ONT R9.4.1 (2026-10-15)

The second part of the benchmark, after the [Col-0 pipeline record](../2026-10-14_benchmark/SUMMARY.md). That
record had six lines of one ecotype, with R10.4.1 reads and no SNP between them. This one has 24 natural
accessions (PRJEB62038) with older R9.4.1 reads: about 88% read identity, base qualities mostly Q3–Q10. They
differ by 9–113 SNPs in the plastome and 1–66 in the mitogenome.

The Illumina truth and BACoN's runs are described in the [minimap2 record](../2026-10-15_minimap2/SUMMARY.md),
which also holds their scripts. Truth: each accession's base at each position, where at least 90% of its
Illumina reads agree at a depth of at least 20; elsewhere `N`. Distances are counted at the positions called in
both accessions.

## Pipelines

| Pipeline | Version | Input | What it gives |
|---|---|---|---|
| BACoN, `--snp-method ska` (default) | 0.3.9 (4967f22) | baited reads; templated or Flye 2.9.6; mitogenome with `--filter-by-reference` | assemblies, distances, VCF |
| BACoN, `--snp-method minimap2` | 96b1808 | the same assemblies | distances, VCF |
| BACoN assemblies, Parsnp | Parsnp 2.1.6 `-c`, as `--snp-method parsnp` | templated plastomes | distances |
| Clair3 | 2.0.3, model r941_prom_hac_g360+g422, `--haploid_precise` | reads aligned with minimap2 to TAIR10.1 (nuclear copy of the mitogenome masked, as for the truth) | SNPs; consensus = reference + SNPs, `N` below 10x |
| [ptGAUL](https://github.com/Bean061/ptgaul) | 1.0.5 (Flye 2.7) | reads; plastome only | assemblies |
| ptGAUL plastomes, MAFFT | 7.526 `--auto`, snp-dists 0.8.2 | ptGAUL plastomes, oriented and rotated to the reference (`orient.py`) | distances |
| ptGAUL plastomes, Snippy | 4.6.0 `--ctgs`, snippy-core | as above | SNPs, distances |
| ptGAUL plastomes, SKA2 | as BACoN compares | as above | distances |
| Snippy, mitogenome | 4.6.0 `--ctgs` | BACoN's Flye mitogenomes | SNPs, distances |

Plastome assembly truth: GetOrganelle 1.7.7.1 on the Illumina reads (`-F embplant_pt -R 15 -k 21,45,65,85,105`).
It gave 22 circular plastomes; Ct-1 and Ge-0 are scaffolds. ptGAUL wrote no path for five accessions where its
graph had a single edge (Are-1, Are-6, Col-0, Ita-0, Jea); its `final_assembly.fasta` was used for those.

## Distances

Against the truth, all 276 pairs ([`scores_distances.tsv`](scores_distances.tsv); matrices in
[`distances/`](distances/)):

| Pipeline | Organelle | Exact | Within 5 | Mean abs. error | Mean signed error | Max error |
|---|---|---|---|---|---|---|
| **Clair3** | plastome | 27 | 210 | **3.8** | −1.2 | 15 |
| BACoN minimap2, Flye | plastome | 12 | 156 | 5.7 | −0.6 | 19 |
| BACoN minimap2, templated | plastome | 6 | 54 | 12.4 | −12.2 | 32 |
| BACoN Parsnp, templated | plastome | 1 | 7 | 33.3 | −33.3 | 56 |
| ptGAUL + SKA2 | plastome | 0 | 10 | 44.7 | −44.6 | 68 |
| BACoN SKA2, Flye | plastome | 1 | 5 | 45.3 | −45.3 | 70 |
| BACoN SKA2, templated | plastome | 1 | 2 | 46.6 | −46.6 | 73 |
| ptGAUL + Snippy | plastome | 0 | 0 | 46.8 | +46.8 | 94 |
| ptGAUL + MAFFT | plastome | 1 | 7 | 902 | +902 | 6,632 |
| BACoN SKA2, Flye | mitogenome | 11 | 97 | **9.6** | −8.5 | 27 |
| BACoN minimap2, templated | mitogenome | 16 | 155 | 10.9 | +9.3 | 42 |
| BACoN SKA2, templated | mitogenome | 6 | 59 | 11.6 | −11.6 | 27 |
| BACoN Flye + Snippy | mitogenome | 3 | 61 | 13.7 | +12.8 | 47 |
| Clair3 | mitogenome | 0 | 20 | 25.1 | +25.1 | 66 |
| BACoN minimap2, Flye | mitogenome | 2 | 19 | 26.7 | −26.6 | 50 |

- **MAFFT** misaligned two plastomes, Col-0 and Jea, which gave every pair with them about 6,600 differences.
  Both are single-edge cases, where ptGAUL's `final_assembly.fasta` is 156 kb. Across the other 22 accessions,
  MAFFT's distances are 51 SNPs too high on average (max 107).
- **Snippy** counted too many plastome SNPs on ptGAUL's plastomes, whose assemblies carry 4–29 substitution
  errors each (below); the cause of the rest of the excess was not investigated.
- **minimap2 with `--ska-min-freq 0.04`** (any one genome covering a position; each pair compared on the
  positions both cover): plastome 7.4 (templated) and 7.1 (Flye), mitogenome 21.4 and 24.5. Better for templated
  plastomes, worse elsewhere.

**At the positions the truth can judge** ([`scores_distances_judgeable.tsv`](scores_distances_judgeable.tsv)):
pairwise differences only at positions called by the truth in both accessions. This is possible only for
pipelines whose sequences are in reference coordinates.

| | Plastome | Mitogenome |
|---|---|---|
| Clair3 | 4.5 | **1.3** |
| BACoN minimap2, templated | 8.6 | 2.1 |
| BACoN minimap2, Flye | **5.6** | 13.3 |
| Snippy | 63.3 | 34.2 |

Most of the mitogenome error is at positions the truth leaves `N`, where the Illumina reads are mixed:
nuclear copies, repeat copies that differ. Clair3 and the templated consensus call a base there, and the
distances count it. At the positions the truth can judge, Clair3 and BACoN's templated mitogenomes with minimap2
are close to the truth (1.3 and 2.1).

## SNPs per accession

Recall (true SNPs found) and precision (calls that are true), at positions the truth calls in that accession
([`scores_snps.tsv`](scores_snps.tsv)):

| Pipeline | Plastome recall / precision | Mitogenome recall / precision |
|---|---|---|
| Clair3 | 0.945 / 0.974 | **0.965 / 0.989** |
| BACoN minimap2 VCF, templated | 0.871 / 0.923 | 0.949 / 0.905 |
| BACoN minimap2 VCF, Flye | 0.913 / 0.851 | 0.800 / 0.699 |
| BACoN SKA2 VCF, templated | 0.741 / **0.990** | 0.704 / 0.988 |
| BACoN SKA2 VCF, Flye | 0.766 / 0.974 | 0.676 / 0.821 |
| Snippy | 0.738 / 0.881 | 0.701 / 0.804 |

- **SKA2:** what it calls is right, but it misses a quarter of the SNPs (the indels next to them).
- **minimap2, templated plastomes:** its 132 false calls fall on only 30 positions, mostly shared by many
  accessions (one by 19). They are local misalignments next to homopolymer deletions of the consensus: minimap2
  places the gap a few bases off and reads the bases beside it as substitutions. Because they are shared, they
  change the distances between accessions little. Masking around indels removed them but lost more true SNPs
  (minimap2 record).

## Plastome assemblies

Against the GetOrganelle truths, with nucmer (`delta-filter -1`); totals over the 24 accessions
([`scores_assembly_*.tsv`](.)):

| Assemblies | Truth covered | Substitutions (max per accession) | Indel events per accession |
|---|---|---|---|
| BACoN templated | all 24 complete | **15 (2)** | 400–488 |
| BACoN Flye (2 contigs each) | 0.997–1 | 69 (8) | 496–687 |
| ptGAUL | 0.998–1 | 275 (29) | 507–742 |
| Clair3 consensus (reference + SNPs) | all 24 complete | 56 (6) | 0–60 (the true indels, which it does not call) |

The templated consensus has the fewest substitution errors; every Nanopore assembly keeps 400–740 homopolymer
indels with these reads. Clair3's consensus has none of those, being the reference with substitutions, but it
misses every true indel.

## Mixed sites in the templated consensus (tested, not kept)

The truth leaves `N` where the Illumina reads are mixed, and the templated mitogenomes call a base there. Two ways
of marking these sites `N` in the consensus were tested:

- **samtools consensus `-A`** (ambiguity codes at heterozygous calls, then `N` as BACoN does; `ambig_run.sh`): at
  most one more `N` per sample. With this error profile, samtools almost never calls a heterozygous site.
- **Masking where the Nanopore reads are mixed:** `N` wherever the major base is below a fraction of the reads,
  on the reference-coordinate consensus (`mask_refcoord.sh`, `mask_ont_counts.py`, `mask_masked.py`). Core
  distances, mean absolute error (signed):

| Major base below | Plastome | Mitogenome | Positions masked per sample (plastome, mitogenome) |
|---|---|---|---|
| (not masked) | 14.3 (−14.3) | 15.0 (+13.8) | — |
| 0.6 | 15.0 (−14.9) | 13.9 (+12.2) | 9, 35 |
| 0.7 | 16.3 (−16.3) | 11.8 (+9.9) | 54, 234 |
| 0.8 | 20.1 (−20.0) | 8.3 (+4.8) | 225, 1,083 |

The mitogenome gains only at 0.8, while the plastome loses true SNPs: with 12% read error, many positions look
mixed. Not kept.

## Conclusions

- **With R9.4.1 reads, SKA2 is the wrong comparison.** The consensus keeps 400–490 indels per plastome, and SKA2
  misses a quarter of each accession's SNPs (recall 0.74) and about half of the SNP sites (242, against 459 with
  minimap2). `--snp-method minimap2` (added after this benchmark) recovers most of
  them: plastome errors of 12.4 (templated) and 5.7 (Flye), against Clair3's 3.8.
- **Mitogenomes:** BACoN's SKA2 comparison of the Flye assemblies has the smallest error (9.6). No pipeline does
  much better, because mixed positions dominate the error. At positions the truth can judge, Clair3 (1.3) and
  BACoN templated with minimap2 (2.1) are both close.
- **Assembly accuracy:** BACoN's templated plastomes have the fewest substitution errors of the assemblers
  (15 in 24 plastomes, against 69 for Flye and 275 for ptGAUL).
- **The other comparison routes are worse:** MAFFT (misalignments and excess SNPs), Snippy (excess SNPs) and Parsnp
  (−33) all did worse than BACoN with minimap2.

## Reproducing

`BENCH` is the benchmark folder. Steps:

1. Truth and BACoN runs: scripts of the [minimap2 record](../2026-10-15_minimap2/SUMMARY.md); its `mm2.py`
   writes the minimap2 comparisons to `acc/compare/minimap2/<run>/`, with `--ska-min-freq 0.04` for the pan
   variant.
2. Plastome truth: `getorganelle.sh`.
3. Clair3: `clair3_run.sh`, then `clair3_distances.sh`.
4. ptGAUL: `ptgaul_run.sh`.
5. Comparisons of assemblies: `assembly_pipelines.sh` (orientation, MAFFT, Snippy, SKA2) and `parsnp.py` on the
   templated plastomes.
6. Scores: `score.py` (distances, judgeable distances, SNPs) and `assembly_score_all.sh` (`assembly_score.py`).
7. Mixed-site test: `ambig_run.sh`, and `mask_refcoord.sh`, `mask_ont_counts.py`, `mask_masked.py`.

Timings are not reported: the pipelines ran at the same time on one machine, with different numbers of samples in
parallel.
