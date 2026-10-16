# `--snp-method minimap2`: SNPs read along the reference (2026-10-15)

A new comparison method, kept as an option (SKA2 stays the default). Each assembly is aligned to the reference
(`minimap2 -x asm20 -c --cs=long`), and its base at every reference position is read from the alignment. Unlike
SKA2, it needs no exact context around a site, so the indels of a consensus do not hide the SNPs next to them.

It was prompted by 24 *Arabidopsis thaliana* accessions sequenced with R9.4.1 reads, with an Illumina truth. There,
SKA2 found plastome distances 46 SNPs too low on average. The cause: the templated plastomes carry 400–490
homopolymer indels each, and a SNP drops out of SKA2's core when an indel lies within 15 bases of it in any of the
24 genomes. R10.4.1 consensus has 15–23 such indels per plastome (Col-0, [mitogenome record](../2026-10-13_mito/SUMMARY.md)),
so it did not show before.

Code: `bacon/compare.py` (`project`, `mask_uncertain`, `run_minimap2`, `minimap2_vcf`). BACoN 0.3.9 with this
change (editable install, environment versions in [`versions.txt`](versions.txt)). The accession assemblies are
from BACoN at commit 4967f22; the Col-0 and potato assemblies are those of the earlier records. Only the
comparison was run again on them ([`scripts/mm2.py`](scripts/mm2.py)).

## How a genome's base is read

At each reference position, a genome has a base, a gap (`-`, a deletion in the genome) or `N`. It is `N` where:

- no alignment covers the position;
- two alignments cover it with different bases (copies of a repeat that differ, or a contig of another origin,
  such as a nuclear copy of the organelle genome, over the organelle's own);
- **at least 3 `N` lie within 50 positions of it** (`mask_uncertain`; the genome's own `N`, or conflicts, but not
  uncovered positions). This rule was added after the first Col-0 results, below.

Insertions in the genomes are left out. `--ska-min-freq` keeps the positions where at least that share of the
sequences (genomes and reference) has a nucleotide; 1 (default) keeps the positions every genome covers. A SNP in
an inverted repeat is counted at both of its positions, like the Illumina truth; SKA2 counts it once.

## 24 *Arabidopsis* accessions, R9.4.1 (PRJEB62038)

**Data and truth.** ONT R9.4.1 whole-genome reads (about 88% read identity; base qualities mostly Q3–Q10) and
Illumina reads of the same 24 accessions. The truth comes from the Illumina reads:

- **Alignment:** the first 10 million read pairs were aligned with bwa mem to TAIR10.1. The nuclear copy of the
  mitogenome on chromosome 2 ([`nuc_mask.bed`](nuc_mask.bed), 259 kb, 95–99% identical) was masked, because its
  reads could not be assigned ([`scripts/acc_truth_illumina.sh`](scripts/acc_truth_illumina.sh)).
- **Base counts:** at every organelle position, at any mapping quality (reads in the inverted repeat and in the
  mitochondrial repeats have MAPQ 0) and base quality ≥ 20 ([`scripts/acc_pileup_calls.py`](scripts/acc_pileup_calls.py)).
  bcftools' genotype calls were not used: at MAPQ 0 they emitted no allele even where every read agreed.
- **Calls:** a sample's base at a position is the one carried by ≥ 90% of its reads at depth ≥ 20, else `N`
  ([`scripts/acc_truth.py`](scripts/acc_truth.py)).

Per accession, the truth calls 154,318–154,469 plastome positions (76 SNPs to the reference on average) and
342,766–367,652 mitogenome positions (57 SNPs). Distances between accessions are 9–113 (plastome) and 1–66
(mitogenome). Tables: [`acc_truth_*.tsv`](.). Distances are compared at the positions called in both accessions.

BACoN ran on the plastome (templated and Flye) and the mitogenome (templated and Flye, `--filter-by-reference`;
[`scripts/acc_bacon.sh`](scripts/acc_bacon.sh)). Scores for all 276 pairs ([`acc_distances.tsv`](acc_distances.tsv)):

| Run | Method | SNP sites | Pairs exact | Within 5 | Mean abs. error | Mean signed error | Max error |
|---|---|---|---|---|---|---|---|
| plastome, templated | SKA2 | 242 | 1 | 2 | 46.6 | −46.6 | 73 |
| plastome, templated | **minimap2** | 459 | 6 | 54 | **12.4** | −12.2 | 32 |
| plastome, Flye | SKA2 | 244 | 1 | 5 | 45.3 | −45.3 | 70 |
| plastome, Flye | **minimap2** | 525 | 12 | 156 | **5.7** | −0.6 | 19 |
| mitogenome, templated | SKA2 | 145 | 6 | 59 | 11.6 | −11.6 | 27 |
| mitogenome, templated | minimap2 | 234 | 16 | 155 | 10.9 | +9.3 | 42 |
| mitogenome, Flye | SKA2 | 174 | 11 | 97 | **9.6** | −8.5 | 27 |
| mitogenome, Flye | minimap2 | 70 | 2 | 19 | 26.7 | −26.6 | 50 |

For scale, Clair3 (R9.4.1 model; haploid calls on the reads aligned to the reference) scored a mean absolute error
of 3.8 on the plastome and 25.1 on the mitogenome. These pipeline comparisons will get their own record.

- **Plastome: a large gain.** minimap2 finds about twice as many SNP sites.
  - Templated assemblies: the remaining shortfall (−12) is SNPs that the templated consensus itself missed.
  - Flye: almost no bias.
- **Mitogenome: no gain.**
  - Templated: the errors are about as large as SKA2's, but now too high. These are sites where the reads are
    mixed: nuclear copies, repeat copies that differ. The consensus picks one allele there, which the truth leaves
    `N`.
  - Flye: the core shrinks to 70 sites. The Flye mitogenomes are fragmented and partly incomplete, and a position
    must be covered by all 24.

## Masking tested, and the rule kept

Tested on the 24 accessions ([`acc_masking_tests.tsv`](acc_masking_tests.tsv);
[`scripts/window.py`](scripts/window.py), [`scripts/window2.py`](scripts/window2.py)):

- **Masking a base whenever an `N` is near it** (`w3`–`w15`): the mitogenome improved little (templated 10.9 to
  10.5), and the plastome lost true SNPs (12.4 to 13.5–14.5).
- **Masking around indels too** (`w*i`): the templated mitogenome improved (6.8 at 15 bases), but the plastome
  fell to 16.5–33.2. R9.4.1 indels are everywhere, true SNPs included.
- **Masking where indels crowd** (`e<W>k<K>`, K indel events within W positions; tried for the potato Flye
  plastomes, below): every setting that changed potato cost true SNPs on the accessions (plastome, templated:
  12.4 to 12.9–22.4).
- **Masking where `N` crowd** (`d<W>k<K>`): at 3 `N` within 20 or 50 positions, the accessions' scores did not
  change (mitogenome, templated: 10.85 to 10.84 at 50). Kept, at 50 positions, for Col-0:

The first Col-0 run of the method found 19 SNPs between the lines in the templated mitogenomes with
`--filter-by-reference`, where the HiFi truth and SKA2 have none. All 19 lie in one line (MA32-119), in 500 bp at
265.3–265.8 kb. There, in every line, the reads are about half and half from two diverged sequences
(copies of a repeat, or a nuclear copy), at the line's usual depth. The consensus calls `N` at most positions; the
few bases it did call followed the other sequence in MA32-119 and the reference in the other lines. The number of false sites by rule: 19 unmasked, 7 at `N` within 3, 3 within 5. With the rule
kept (3 `N` within 50): none.

## Col-0 lines, R10.4.1

Six mutation-accumulation lines of one ecotype ([`col0_lines.tsv`](col0_lines.tsv)):

| Run | SKA2 sites / max between lines | minimap2 sites / max between lines |
|---|---|---|
| plastome, all four runs | 0 / 0 | 0 / 0 |
| mitogenome, templated, `--filter-by-reference` | 0 / 0 | 0 / 0 |
| mitogenome, Flye, `--filter-by-reference` | 1 / 1 | 0 / 0 |
| mitogenome, templated, default filter | 38 / 8 | 66 / 8 |
| mitogenome, Flye, default filter (the plastome was assembled) | 28 / 0 | 49 / 0 |

The default-filter mitogenomes were already wrong (the [mitogenome record](../2026-10-13_mito/SUMMARY.md)); with
minimap2 their distances to the reference are larger.

## Simulated scenarios

`run_validation.sh` now also compares the templated and Flye assemblies with minimap2 ([`SCORES_*.md`](.)).
Exact SNP distances:

| Scenario | Templated, SKA2 | Templated, minimap2 | Flye, SKA2 | Flye, minimap2 |
|---|---|---|---|---|
| plastid | 45/45 | 45/45 | 45/45 | 45/45 |
| linear | 10/10 | 10/10 | 10/10 | 10/10 |
| rdna | 10/10 | 10/10 | 1/10 | **6/10** |
| heteroplasmy | 45/45 | 36/45 | 28/45 | 28/45 |
| numt | 11/21 | 11/21 | **16/21** | 10/21 |

minimap2's best convention in the scores is `sites` for plastid and rdna: it counts each changed position, so a
SNP in both copies of the inverted repeat counts twice.

- **heteroplasmy, templated:** at an even mixture (h050), minimap2 counts one SNP (truth 0). The consensus called
  the minor allele in one copy of the inverted repeat (position 44,351) and not in the other, and SKA2 does not
  count such a site.
- **heteroplasmy, Flye:** at 0.4 and 0.5, minimap2 counts 5 and 10 SNPs, minor alleles that the assemblies carry;
  SKA2 counts 2 and 5.
- **numt, Flye:** contigs made from nuclear copies of the plastome overlap the plastome contig in n04 and n06. Where
  their bases differ, the position becomes `N`, which takes three true SNPs out of the core (one in n04, two in
  n06). Two consensus errors of n04's Flye plastome count as SNPs; SKA2 skipped them.

Masking: the rule kept changed none of these distances (each run compared with and without it).

## Potato, R9.4.1 (28 cultivars)

No truth; the agreement of templated and Flye assemblies of the same cultivars, 406 pairs with the reference
([`potato_templated_vs_flye.tsv`](potato_templated_vs_flye.tsv)). The assemblies are those of the
`--filter-by-reference` plastome runs and of the mitogenome runs of the
[mitogenome record](../2026-10-13_mito/SUMMARY.md).

| | Identical pairs | Largest difference | Mean difference |
|---|---|---|---|
| plastome, SKA2 | 326 | 2 | 0.2 |
| plastome, minimap2 | 48 | 195 | 28.9 |
| mitogenome, SKA2 | 379 | 6 | 0.1 |
| mitogenome, minimap2 | 77 | 162 | 13.4 |

The Flye plastomes disagree because minimap2 counts every difference they carry. In Bankir, a single circular
contig, the last 8.6 kb of the alignment (reference 146.7–155.3 kb) has 149 indels and 142 substitutions (about
1.7 of each per 100 bases), against 124 indels and 2 substitutions in the other 147 kb. This unpolished stretch
gives Bankir 113 more SNPs to the reference than its templated assembly. SKA2 finds no split k-mer in such a
stretch and skips it. Masking crowded indels did not close the gap without costing true SNPs on the accessions
(above).

## Decision

`--snp-method minimap2` is kept as an option, recommended for templated assemblies of older Nanopore reads
(R9.4.1). On those, it cut the plastome distance error from 46.6 to 12.4 SNPs, and it matched SKA2 on the
simulated scenarios and the Col-0 lines. SKA2 stays the default:

- minimap2 is no better on mitogenomes;
- with de novo assemblies, it counts their uncorrected errors (potato);
- on R10.4.1 data, SKA2 was already exact.

## Erratum: the Col-0 pipeline record

The [Col-0 pipeline record](../2026-10-14_benchmark/SUMMARY.md) gives BACoN's version as commit 7c98957. Its
timing runs actually used a BACoN installed on 2026-10-13 at 19:46 from an uncommitted working tree. That tree
still held the competitive-baiting code, which was tested and dropped (it adds a `Decoy_reads` column) and is
inactive unless its option is given. The read selection, the assemblies and the distances of those runs were
checked to be the same as 7c98957's. The environment has since been reinstalled from 4967f22.

## Reproducing

- **Accession truth:** `acc_truth_illumina.sh`, then `acc_pileup.sh` (base counts and truth tables).
- **BACoN runs:** `acc_bacon.sh`.
- **minimap2 on existing runs:** `mm2.py OUT RUN_DIR`.
- **Masking tests:** `window.py` and `window2.py`.
- **This record's tables:** `tables.py`.
- **Simulations:** `run_validation.sh`, then `rerun_mm2.sh`. The first pass ran with an earlier version of the
  method; `rerun_mm2.sh` runs the minimap2 comparisons again with the final code.

`BENCH` is the benchmark folder: `acc/`, `col0/`, `potato_fref/`, `potato_mito/`, `refs/`.
