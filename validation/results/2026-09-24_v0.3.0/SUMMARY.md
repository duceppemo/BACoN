# Validation of BACoN 0.3.0 (2026-09-24)

Two parts: a **screening** of candidate programs, which decided the programs BACoN 0.3.0 uses, and the
**validation** of the released code with the chosen programs. Simulated data only in this folder; the
real-data comparison (public potato data) is summarized at the end.

## Validation of the released code

```bash
bash validation/run_validation.sh <folder> 16      # BACoN 0.3.0, conda environment from environment.yml
```

Program versions: [`versions.txt`](versions.txt) (minimap2 2.31, Filtlong 0.3.1, Flye 2.9.6, myloasm 0.7.0,
samtools 1.24, SKA2 0.5.1, Parsnp 2.1.5, FastTree 2.2.0). Scores: [`SCORES_plastid.md`](SCORES_plastid.md),
[`SCORES_linear.md`](SCORES_linear.md), [`SCORES_rdna.md`](SCORES_rdna.md) (and `.json`); wall times:
[`timings.txt`](timings.txt). Scenarios: [`../../README.md`](../../README.md).

### Assemblies (errors = mismatches + inserted + deleted bases against the true genome)

| Scenario | Assembly | Samples | Contigs | Circular | Errors | Notes |
|---|---|---|---|---|---|---|
| plastid (9 samples) | samtools (default) | 9/9 | 9 | NA | 1 | the one mismatch is in the 10x sample |
| | Flye | 9/9 | 11 | 8 | 68 | all 68 in the 10x sample (3 contigs, 19 kb duplicated); 0 in the others |
| | myloasm | 9/9 | 10 | 8 | 0–1 | 10x sample: 2 contigs, 9 kb duplicated |
| linear (4 samples, 3 molecules each) | samtools | 4/4 | 12 | NA | 0 | |
| | Flye | 4/4 | 12 | 0 | 0 | no linear molecule circularized; ends 4% short |
| | myloasm | 4/4 | 12 | 0 | 2 | |
| rdna (4 samples) | samtools | 4/4 | 4 | NA | 0 | |
| | Flye | 4/4 | 4 | **4** | 0 | each array collapsed into one 6 kb **circular** contig (27% of the molecule; flanks lost) |
| | myloasm | 4/4 | 4 | 0 | 3001 | 3 of 4 correct; one sample lost one 3 kb copy of the unit |

### SNP distances (pairs of genomes with exactly the true distance)

| Scenario | samtools + SKA2 (default) | Flye + SKA2 | myloasm + SKA2 | Flye + SKA2 pan (0.5) | Flye + Parsnp `-c` | Flye + Parsnp default |
|---|---|---|---|---|---|---|
| plastid | **45/45** | 45/45 | 45/45 | 45/45 | 36/45 (misses SNPs in the inverted repeat; off by 1 next to the deletion) | 22/36, **drops the 10x sample** |
| linear | **10/10** | 10/10 | 6/10 (max error 2) | 10/10 | 10/10 | 10/10 |
| rdna | **10/10** | 1/10 (flanks lost) | 6/10 (max error 1) | 1/10 | fails: Parsnp rejects assemblies 3.7 times shorter than the reference | — |

SKA2 counts a SNP in an inverted repeat or in every copy of a tandem unit once (the `events` convention);
Parsnp misses SNPs in inverted repeats (`no_repeat`).

Wall time, plastid scenario (9 samples, `-t 16 -p 4`): samtools 6.9 s, myloasm 11.8 s, Flye 60.3 s.

**myloasm is not fully deterministic**: in two runs on byte-identical reads, the assembly of the 10x plastid
sample differed by one inserted base (all other results identical).

## Screening of candidate programs

Run the same day with a development version of BACoN that also supported Shasta 0.14.0, Raven 1.8.3,
Rebaler 0.2.0, Medaka 2.2.2 (polishing, and templated assembly with gaps filled with `N` or with the
reference) and Snippy 4.6.0; those were removed afterwards. Scores in [`screening/`](screening/) (the rdna
scenario on four assemblers only; the myloasm scores there predate the SKA2 fix below). Main results on the
plastid scenario:

| Candidate | Consensus errors (9 samples) | SNP distances exact (SKA2) | Time | Notes |
|---|---|---|---|---|
| samtools consensus (templated) | 1 | 45/45 | 4 s | kept, default |
| Flye | 68 (10x sample) | 45/45 | 11 s | kept |
| myloasm | 0 | 45/45 after the fix below (21/45 before) | 5 s | kept |
| Shasta | 20 | 36/45 | 5 s | ~20 bp unclipped overlaps at circular ends; lost 3.6 kb of an inverted repeat in one sample; fragmented the rdna arrays; on real data, some assemblies were a fraction of the genome and no core SNP remained. Removed |
| Raven | 138 | 3/45 (max error 35) | 5 s | homopolymer deletions; **no contig at all** for linear molecules under 10 kb. Not adopted |
| Rebaler (templated, BACoN 0.2) | 57 | 36/45 (max error 24) | 20 s | unmaintained, needs `setuptools<81`. Replaced by samtools consensus |
| Medaka templated | 38 | 45/45 | 55 s | worse consensus than samtools, 8x slower. Not adopted |
| Flye + Medaka polishing | 71 | 45/45 | 86 s | no gain on simulated data; on real potato data it moved 80 of 406 distances away from the value on which the templated and the Flye assemblies agree. Not adopted |
| Snippy (on Flye) | — | 45/45 only without the IR SNPs | — | misses SNPs in inverted repeats; cannot be installed with current Flye and Filtlong. Removed |
| Parsnp default (on Flye) | — | dropped the 10x sample | — | BACoN uses `-c` |

On the linear scenario Flye circularized nothing; on the rdna scenario it collapsed every tandem array into a
circular unit (Shasta fragmented them; myloasm kept them linear).

## Bugs found by the validation (all fixed before the validation of the released code)

- myloasm's circular contigs were not recognized (`circular-yes_` in its headers), and circular contigs lost
  the SNPs within 30 bp of their start in SKA2, whose split k-mers do not wrap around: circular contigs are
  now extended by 30 bases before SKA2 (myloasm plastid: 21/45 → 45/45 exact distances).
- When no SNP site is shared by all genomes, FastTree failed and the run stopped without writing
  `summary.tsv`: BACoN now skips the tree with an explanation, and always writes the summary.
- Reads given as fasta failed at filtering (Filtlong needs qualities).
- `--template-gaps reference` filled an all-`N` consensus twice.
- Scoring: intervals of different molecules were merged, and repeated regions were counted as duplicated
  (fixed in `score.py` before any score here was kept).

## Real data (public)

The 28 potato cultivars of BioProject PRJNA807056 (MinION, R9.4.1-era reads, mean Q14–16; reference
NC_008096.2), analysed as in the [tutorial](https://github.com/duceppemo/BACoN/wiki/Tutorial), with
`validation/compare_runs.py`:

| Comparison | Pairs with identical SNP distances (of 406) | Largest difference |
|---|---|---|
| samtools + SKA2 vs Flye + SKA2 | **406** | 0 |
| Flye + SKA2 vs Flye + Medaka + SKA2 | 326 | 4 |
| Flye + SKA2 vs Flye + Parsnp `-c` | 121 | 10 |
| Flye + SKA2 vs Flye + Snippy | 97 | 22 |
| Flye + SKA2 vs myloasm + SKA2 | 93 | 21 |

The templated and the de novo assemblies, analysed independently, agree exactly. Flye assembled 17 of the
28 plastomes as one circular contig (myloasm: 8); Flye resolved the 241 bp *ndhC*–*trnV* insertion that
separates the cytoplasm types exactly, while the templated consensus placed it with some bases called `N`.
