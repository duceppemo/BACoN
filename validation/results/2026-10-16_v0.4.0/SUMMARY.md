# BACoN 0.4.0 on the simulated scenarios, and a stricter consensus cutoff (2026-10-16)

Two parts:
1. The simulated validation rerun with the released 0.4.0 (tag v0.4.0, installed from its tarball), whose default
   comparison is now minimap2 for templated assemblies and SKA2 for de novo ones.
2. A change after 0.4.0: the templated consensus is called with `samtools consensus -C 20` instead of samtools'
   default 10. A base called with a quality below 20 becomes `N`.

## 1. Simulated scenarios with 0.4.0

`run_validation.sh` (now: each assembler with its default comparison, then the other method on the same
assemblies). Scores: [`SCORES_*.md`](.); program versions: [`versions.txt`](versions.txt). Exact SNP distances:

| Scenario | Templated, minimap2 (default) | Templated, SKA2 | Flye, SKA2 (default) | Flye, minimap2 | myloasm, SKA2 | Flye, Parsnp `-c` |
|---|---|---|---|---|---|---|
| plastid | 45/45 | 45/45 | 45/45 | 45/45 | 45/45 | 36/45 |
| linear | 10/10 | 10/10 | 10/10 | 10/10 | 6/10 | 10/10 |
| rdna | 10/10 | 10/10 | 1/10 | 6/10 | 6/10 | fails |
| heteroplasmy | 36/45 | 45/45 | 28/45 | 28/45 | — | — |
| numt | 11/21 | 11/21 | 16/21 | 10/21 | — | — |

Every value equals the one obtained with the development versions in the
[minimap2 record](../2026-10-15_minimap2/SUMMARY.md) and the 0.3.8 record. Parsnp still fails on rdna, as before.

## 2. Consensus cutoff 20 (after 0.4.0)

`samtools consensus` (Bayesian mode, `-X r10.4_sup`) gives each call a quality, and calls `N` below the cutoff,
10 by default. With few reads, a base can be called from a mixture of bases at a quality just above 10, and the
minimap2 comparison then counts it as a SNP, which SKA2 would often have skipped with its neighbours.

**Real data.** [`scripts/consc.py`](scripts/consc.py) rebuilt each templated consensus from the alignment BACoN kept
(`3_assembled/<sample>/aligned.bam`), with BACoN's options and `-C` 10, 15, 20 or 30. It then compared the
assemblies with minimap2, as `--snp-method minimap2` does. Cutoff 10 reproduced the published SNP sites exactly (459,
234, 447, 0, 0, 170). Scores ([`consensus_cutoff_real.tsv`](consensus_cutoff_real.tsv),
[`scripts/cutoff_score.py`](scripts/cutoff_score.py)):

| Data (truth) | `-C 10` | `-C 15` | `-C 20` | `-C 30` |
|---|---|---|---|---|
| *Arabidopsis* plastomes, R9.4.1 (Illumina): mean abs. error, SNP sites | 12.40, 459 | 12.40, 459 | 12.76, 458 | 12.76, 458 |
| *Arabidopsis* mitogenomes, R9.4.1, `--filter-by-reference` (Illumina) | 10.85, 234 | 10.47, 232 | 10.41, 230 | 10.29, 229 |
| human mitogenomes, R10.4.1 (Illumina) | 0.87, 447 | 0.87 | 0.87 | 0.87 |
| Col-0 plastomes and mitogenomes, R10.4.1 (0 between lines) | 0 | 0 | 0 | 0 |
| potato plastomes, R9.4.1 (no truth): distinct plastomes | 7 | 7 | 7 | 7 (Otrada 1 SNP from its group, not 2) |

**Simulations.** The full validation with `-C 20` ([`consensus_cutoff_20/`](consensus_cutoff_20/)) gave the same
templated SNP distances as 0.4.0 in every scenario. Templated assemblies had more `N` where support is weak:
- the 10x plastid sample: 23 `N` (1 with 0.4.0), so 21 mismatches against the truth, `N` included (1);
- heteroplasmy: 6 `N` at 0.3 (5 with 0.4.0), 23 at 0.4 (22);
- numt: errors inside the copied segments, nearly all `N`: 1,276 in n04 (1,230), 25 in n06 (22).

The Flye rows differ too, from Flye's own run-to-run variation (the cutoff applies only to the templated
consensus).

**Decision: kept, at 20.** A base with little support becomes `N` rather than a SNP. The *Arabidopsis* mitogenome
distances improved (10.85 to 10.41). The plastome distances lost one SNP site (459 to 458; error 12.40 to 12.76).
The human and Col-0 results and every simulated distance were unchanged. 30 was not chosen: it gained little over
20 on the mitogenomes and lost the same plastome site. The cutoff is part of the templated assembly's checkpoint: a
folder of an earlier version remakes its templated assemblies on the next run (baiting and filtering are kept).
