# Human mitogenomes, ONT R10.4.1: 32 samples of the 1000 Genomes long-read cohort (2026-10-16)

The third part of the benchmark ([Col-0 record](../2026-10-14_benchmark/SUMMARY.md),
[accession record](../2026-10-15_accessions/SUMMARY.md)): many genetically different mitogenomes of one species,
with current reads and a truth from the same cell lines.

## Data

- **Nanopore:** 32 samples of the 1000 Genomes Long-Read Sequencing Consortium (S3 `1000g-ont`), all SQK-LSK114
  R10.4.1 ([`samples_ont.tsv`](samples_ont.tsv)).
  - Ancestry: 6 per continental group (8 East Asian), including a mother and child, HG02025 and HG02024
    (pedigree of the 1000 Genomes project).
  - Basecallers: Dorado SUP v5.0.0 for 23 samples, v4.3 for one, v3.3 for three; not stated in the file name for
    five.
  - Reads: those aligned to chrM in each sample's public hg38 alignment, primary alignments only
    ([`scripts/fetch_chrM.sh`](scripts/fetch_chrM.sh)). 2,388–106,310 reads per sample, 976–6,251x on chrM.
  - These reads were already selected by alignment to the whole genome, so this tests BACoN's assembly and
    comparison, not its baiting: reads of nuclear copies of mtDNA mostly aligned to the nucleus.
- **Truth:** the 1000 Genomes 30x Illumina alignments of the same cell lines, region chrM (about 14,000x).
  - Calls: base counts at any mapping quality and base quality ≥ 20; a sample's base is the one carried by ≥ 90%
    of its reads at depth ≥ 20, else `N`, as for the accessions ([`scripts/truth.py`](scripts/truth.py)).
  - Coverage: 16,563–16,568 of 16,569 positions called per sample ([`truth_callable.tsv`](truth_callable.tsv)).
  - Variation: 10–94 SNPs to the reference, and 452 variable positions across the 32 samples. Distances between
    samples are 0–105 (mean 46).
  - Pairs that must be identical: the mother and child have a distance of 0 in the truth, as expected. So do
    HG00731 and HG01061 (both Puerto Rican; no relation in the pedigree).
- **Reference:** rCRS, NC_012920.1. It is identical to the chrM the Illumina reads were aligned to (same MD5).
- **The GIAB samples planned for this part were dropped:** the remote fetch of their chrM reads failed partway
  through the BAMs (HG002 to HG007), and none has an Illumina truth here.

## Pipelines

- **BACoN at 4d406e8:**
  - assembly: templated (default) and Flye;
  - comparison: each assembly set compared with SKA2 (default) and with `--snp-method minimap2`;
  - run with `-t 32 -p 8` ([`scripts/bacon_run.sh`](scripts/bacon_run.sh)).
- **Clair3 2.0.3** ([`scripts/clair3_run.sh`](scripts/clair3_run.sh)):
  - input: the reads aligned to rCRS with minimap2;
  - options: `--haploid_precise`, with the R10.4.1 SUP model matching the basecaller (v5.0.0 or not stated:
    `r1041_e82_400bps_sup_v500`; v3.3 and v4.3: `sup_v410`, the closest available);
  - consensus: rCRS with the PASS substitutions, `N` below 10x.

## Distances

Against the truth, all 496 pairs ([`scores_distances.tsv`](scores_distances.tsv)):

| Pipeline | SNP sites | Exact | Within 1 | Within 5 | Mean abs. error | Mean signed error | Max error |
|---|---|---|---|---|---|---|---|
| **Clair3** | — | 320 | 436 | 496 | **0.54** | −0.03 | 4 |
| **BACoN templated, minimap2** | 447 | 193 | 406 | 496 | **0.87** | +0.41 | 4 |
| BACoN Flye, minimap2 | 388 | 3 | 13 | 219 | 6.9 | −6.9 | 23 |
| BACoN templated, SKA2 (default) | 174 | 2 | 2 | 2 | 31.0 | −31.0 | 68 |
| BACoN Flye, SKA2 | 144 | 2 | 2 | 2 | 33.8 | −33.8 | 74 |

Every pipeline gives the mother and child a distance of 0. Clair3 and minimap2 give HG00731 and HG01061 a distance
of 2. HG00731 is heteroplasmic at both positions: 82–86% of its Illumina reads carry the other base, so the truth
leaves them `N`, and the Nanopore pipelines call the major allele.

**At the positions the truth can judge** ([`scores_distances_judgeable.tsv`](scores_distances_judgeable.tsv);
pairwise differences only at positions the truth calls in both samples):

| | Exact | Within 1 | Mean abs. error | Max error |
|---|---|---|---|---|
| BACoN templated, minimap2 | 406 | 464 | **0.25** | 3 |
| Clair3 | 406 | 464 | 0.31 | 4 |
| BACoN Flye, minimap2 | 379 | 420 | 0.53 | 11 |

## SNPs per sample

Recall and precision at the positions the truth calls in each sample ([`scores_snps.tsv`](scores_snps.tsv)):

| Pipeline | True SNPs found (of 1,181) | Recall | Precision |
|---|---|---|---|
| BACoN templated, minimap2 VCF | 1,177 | **0.997** | 1.0 |
| Clair3 | 1,176 | 0.996 | 1.0 |
| BACoN Flye, minimap2 VCF | 1,174 | 0.994 | 1.0 |
| BACoN templated, SKA2 VCF | 1,018 | 0.862 | 1.0 |
| BACoN Flye, SKA2 VCF | 1,016 | 0.860 | 1.0 |

No pipeline made a false call at a position the truth calls.

## Why SKA2 fails on these mitogenomes

The cause is not read errors: the templated consensus differs from the truth at almost no position. It is the
density of variation. Of the 452 variable positions in these 32 genomes, 272 have another variant within 15 bases
in at least one sample; 100 are in the control region (16024–576).

SKA2 keeps a core site only when its 15-base flanks are identical in all 33 sequences. That left 174 of the
452 sites, against 447 for minimap2. Its VCF (from `ska map`, each genome against the reference) needs only that
genome's own flanks to match the reference, and found 86% of the SNPs.

The earlier R10.4.1 data had few, scattered SNPs: the Col-0 lines differ by none, and the simulated genomes by
1–40 spread over 70 kb. That is why the failure did not show before.

## Flye assemblies

Flye assembled 28 of the 32 mitogenomes as one circular contig ([`flye_assemblies.tsv`](flye_assemblies.tsv)).
Eight of them are about twice the reference length (1.94–2.00x): the genome twice over. BACoN's note flags them.

Otherwise:
- NA19429: one contig of 0.83x;
- HG02601: 4 contigs;
- NA19238, HG01063, NA20847: one contig each of 1.15–1.27x, not circular.

minimap2 reads a doubled genome as one when both copies carry the same bases. Its Flye distances are 6.9 SNPs
too low because two assemblies remove 61 sites from the core, for every pair:
- NA19429's incomplete assembly removes 36;
- HG02024's doubled assembly removes 25: its two copies differ there, so those positions become `N`.

At judgeable positions, compared pair by pair, the error is 0.53.

## Time

The full templated run took 59 s for all 32 samples (32 threads, 8 samples at a time); the Flye run took 570 s.
Each comparison took 2 s.

## Conclusions

- **BACoN with the templated assembly and `--snp-method minimap2` matches Clair3 on human mitogenomes:**
  - judgeable distance error 0.25 (Clair3 0.31);
  - SNP recall 0.997 with no false call.
- **SKA2, the default, is wrong here:** distances are 31 SNPs too low on average, although the reads are current
  R10.4.1. Where genomes differ at many nearby positions, SKA2's core loses most of the sites.
- **Datasets where SKA2 was right:** genomes with few, scattered SNPs (the Col-0 lines, the simulations).

## Erratum: the minimap2 record

The [minimap2 record](../2026-10-15_minimap2/SUMMARY.md) kept SKA2 as the default partly because "on R10.4.1 data,
SKA2 was already exact". That holds only for genomes with few, scattered SNPs. On these 32 human mitogenomes
(R10.4.1), SKA2's distances were 31 SNPs too low, and minimap2's were within 1 SNP for 406 of 496 pairs.

## Reproducing

`BENCH` is the benchmark folder. Steps:

1. `fetch_chrM.sh`: the reads and the Illumina alignments; its GIAB part failed (above).
2. `pileup.sh`, then `truth.py`: the truth.
3. `bacon_run.sh`: the BACoN runs.
4. `clair3_run.sh`, then `clair3_distances.sh`: the Clair3 baseline.
5. `score.py`: the scores.
