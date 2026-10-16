# Validation

BACoN 0.3 was rebuilt around one rule: a program is used only if it is maintained and does at least as well
as the alternatives on data where the answer is known. The records, with every score, are in
[`validation/results/`](https://github.com/duceppemo/BACoN/tree/main/validation/results); this page
summarizes the [2026-09-24 record](https://github.com/duceppemo/BACoN/blob/main/validation/results/2026-09-24_v0.3.0/SUMMARY.md)
(BACoN 0.3.0). The [2026-10-13 record](https://github.com/duceppemo/BACoN/blob/main/validation/results/2026-10-13_v0.3.8/SUMMARY.md)
reran it with BACoN 0.3.8 and the same programs: every score on simulated data was identical, and on the potato
data the templated and the Flye assemblies still agree on all 406 SNP distances.

## Mitogenomes, nuclear copies and heteroplasmy

The [mitogenome record](https://github.com/duceppemo/BACoN/blob/main/validation/results/2026-10-13_mito/SUMMARY.md)
adds real mitogenomes with a truth, and two simulated scenarios (`heteroplasmy`, `numt`).

On six *Arabidopsis* Col-0 lines (ONT R10.4.1 whole-genome reads; truth: each line's organelles assembled from its
PacBio HiFi reads):

| | Substitutions vs truth | `N` bases | SNPs to the reference | Max SNPs between lines |
|---|---|---|---|---|
| plastome, templated or Flye, with or without `--filter-by-reference` | 0 | 0 | 0 | 0 |
| mitogenome, templated | 189–289 | 19,227–278,015 | 29–34 | 8 |
| mitogenome, templated, `--filter-by-reference` | **0** | **105–162** | **0** | **0** |
| mitogenome, Flye | the plastome was assembled (94–96% off the reference) | — | — | — |
| mitogenome, Flye, `--filter-by-reference` | 0–6 | — | 0–1 | 1 |

By default, Filtlong filled the target depth with long reads of which only a stretch matched the mitogenome;
`--filter-by-reference` ranks the reads by their k-mers found in the reference ([Methods](Methods#2-filtering)).
It is not the default: on the potato plastomes, one cultivar gained a 1 kb block of `N` and fewer Flye plastomes
were circular (11 instead of 17), with the same Flye distances. `Off_reference_pct` stayed at 0–0.3% for the
simulated assemblies and the potato and Col-0 plastomes, at 0–8.3% for the Col-0 mitogenomes assembled with
`--filter-by-reference`, and flagged every Flye "mitogenome" that was the plastome. The simulated
heteroplasmy and nuclear copies are summarized in [Methods](Methods#heteroplasmy).

## Other organelle pipelines

The [benchmark record](https://github.com/duceppemo/BACoN/blob/main/validation/results/2026-10-14_benchmark/SUMMARY.md)
ran four other pipelines that take Nanopore reads on the same six Col-0 lines, scored against the same truths:

| Pipeline | Plastomes complete, substitutions | Mitogenomes complete, substitutions | Max SNPs between lines (plastome, mitogenome) |
|---|---|---|---|
| BACoN, templated (mitogenome with `--filter-by-reference`) | 6, 0 | 6, 0 | 0, 0 |
| BACoN, Flye (mitogenome with `--filter-by-reference`) | 6, 0 | 6, 0–6 | 0, 1 |
| ptGAUL | 6, 0 | — | 0, — |
| TIPPo | 6, 0 | 4 (Flye failed on two lines), 0–18 | 0, 1 |
| Oatk (k 501) | 4 (two in 139 and 184 fragments), 0 | 0 (11–70% of the truth) | 0, 1 |
| CLAW, Flye | 5 (one 43 kb fragment), 0 | — | 0, — |
| CLAW, Unicycler | 5 (three with extra contigs, 482–682 kb), 0–176 | — | 0, — |

Every pipeline's complete plastomes had 9–51 indels against the truth (in BACoN's, nearly all homopolymer runs
one or two bases short), and 222 in one of CLAW-Unicycler's; none had substitutions but CLAW-Unicycler's. Alone
on the machine with 16 threads, per line: BACoN 9–15 min per organelle and 1.3 GB, ptGAUL 14 min and 7.6 GB, CLAW
30 min and 4.4 GB, Oatk 46 min and 40 GB, TIPPo 59 min and 85 GB (Oatk and TIPPo build both organelles in one
run). The record also corrects the share of the truth covered given in the mitogenome record (an erratum: complete
genomes cover 0.998–1.0 of their truth, not 0.6–0.99).

## SNP methods: SKA2 and minimap2

The [minimap2 record](https://github.com/duceppemo/BACoN/blob/main/validation/results/2026-10-15_minimap2/SUMMARY.md)
adds `--snp-method minimap2`, tested on 24 *Arabidopsis* accessions with ONT R9.4.1 reads (PRJEB62038). The truth
is each accession's Illumina reads: the base carried by at least 90% of the reads at a depth of at least 20. Mean
absolute error of the 276 pairwise distances:

| | SKA2 | minimap2 |
|---|---|---|
| plastome, templated | 46.6 (all too low) | **12.4** |
| plastome, Flye | 45.3 | **5.7** |
| mitogenome, templated, `--filter-by-reference` | 11.6 | 10.9 |
| mitogenome, Flye, `--filter-by-reference` | **9.6** | 26.7 |

- **Why SKA2 fails on the plastomes:** the consensus of R9.4.1 reads keeps 400–490 homopolymer indels per
  plastome, and SKA2 drops a site when an indel lies within 15 bases of it in any genome.
- **Where minimap2 equals SKA2:** on the Col-0 lines (R10.4.1; no SNP between lines) and on the simulated
  plastid, linear and rdna scenarios, where minimap2 was also better on the Flye rdna assemblies.
- **Where minimap2 counts errors SKA2 skips:** on the potato plastomes, the Flye and templated distances agree on
  48 of 406 pairs (SKA2: 326), because minimap2 counts the unpolished stretches of some Flye assemblies.

So minimap2 is an option, recommended for templated assemblies of R9.4.1 reads; SKA2 stays the default.

## Simulated data with a known truth

[`validation/simulate.py`](https://github.com/duceppemo/BACoN/blob/main/validation/simulate.py) generates
three scenarios, with reads that mimic R10.4.1 SUP reads (substitutions, indels, homopolymer errors) and 15%
off-target reads:

- **plastid**: a 70 kb chloroplast-like genome with two inverted repeats, and nine samples: two identical to
  the reference, two 1–2 SNPs away, a clade of three, one with 40 SNPs and a 1.5 kb deletion, one at 10x depth;
  SNPs inside the inverted repeat and a 5 bp insertion.
- **linear**: three linear molecules of 15, 9 and 6 kb.
- **rdna**: a linear molecule with a 3 kb unit repeated six times between two flanks, like an rDNA array.

Results with the programs of BACoN 0.3:

| | samtools (templated, default) | Flye | myloasm |
|---|---|---|---|
| Consensus errors, plastid, 8 samples at 40x | 0 | 0 | 0 |
| Plastid sample at 10x | 1 contig, 1 error | 3 contigs, 68 errors | 2 contigs, 0–1 errors (not deterministic) |
| Linear molecules | correct | correct, none circularized | 2 errors in one sample |
| Tandem array (rdna) | correct | **collapsed into a circular 6 kb contig** | correct in 3 of 4 samples |
| Exact SNP distances with SKA2 (plastid, linear, rdna) | 45/45, 10/10, 10/10 | 45/45, 10/10, 1/10 | 45/45, 6/10, 6/10 |

| SNP method (on the Flye assemblies, plastid) | Exact distances | Remark |
|---|---|---|
| SKA2 core (default) | 45/45 | counts a SNP in the inverted repeat once |
| SKA2 pan (`--ska-min-freq 0.5`) | 45/45 | |
| Parsnp `-c` | 36/45 | misses SNPs in the inverted repeat; off by one next to the deletion |
| Parsnp default | 22/36 | dropped the 10x sample |

## Programs considered and not kept

| Program | Result |
|---|---|
| Shasta | overlapping ends on circular contigs, part of an inverted repeat lost, fragmented tandem arrays; on real data some assemblies were a fraction of the genome |
| Raven | homopolymer deletions; no contig at all for linear molecules under 10 kb |
| Rebaler | 57 consensus errors on the plastid scenario and SNP distances wrong by up to 24; unmaintained |
| Medaka (polishing, or templated) | no gain on simulated data; on real data it moved 80 of 406 distances away from the value on which the templated and de novo assemblies agree; slower, needs PyTorch and a model matching the basecaller |
| Snippy | misses SNPs in inverted repeats; cannot be installed with current assemblers |
| PhaME | its dependency check compares versions as decimals and rejects samtools 1.10 to 1.29 |

## Real data

The 28 potato cultivars of the [Tutorial](Tutorial) (public; older MinION reads, mean Q14–16): the templated assembly and the
Flye assembly, analysed independently, give **identical SNP distances for all 406 pairs** of genomes (28 cultivars and the reference). Both
separate the T-type plastomes from the others, which carry the 241 bp sequence absent from T-type plastomes
(resolved exactly by Flye).

## Reproducing

```bash
conda activate BACoN
bash validation/run_validation.sh /tmp/bacon_validation 16
```

The simulated reads are generated from a fixed seed (byte-identical files). All results reproduce exactly,
except that myloasm is not fully deterministic: its assembly of the 10x sample differed by one base between
two runs.
