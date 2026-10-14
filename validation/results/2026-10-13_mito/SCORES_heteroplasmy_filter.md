## Scenario `heteroplasmy`

### Assemblies

| Run | Sample | Contigs | Circular (truth) | Length / truth | Covered | Duplicated bp | Mismatches | Ins bp | Del bp | Errors /100 kb |
|---|---|---|---|---|---|---|---|---|---|---|
| flye | h000 | 1 | 1 (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | h002 | 1 | 1 (1) | 70,001 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | h005 | 1 | 1 (1) | 69,999 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | h010 | 1 | 1 (1) | 69,999 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | h020 | 1 | 1 (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | h030 | 1 | 1 (1) | 69,999 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | h040 | 1 | 1 (1) | 70,000 / 70,000 | 1.0000 | 0 | 3 | 0 | 0 | 4.3 |
| flye | h050 | 1 | 1 (1) | 70,001 / 70,000 | 1.0000 | 0 | 8 | 1 | 0 | 12.9 |
| flye | h020s | 1 | 1 (1) | 70,001 / 70,000 | 1.0000 | 1 | 0 | 0 | 0 | 0.0 |
| samtools | h000 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | h002 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | h005 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | h010 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | h020 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | h030 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 5 | 0 | 0 | 7.1 |
| samtools | h040 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 21 | 0 | 0 | 30.0 |
| samtools | h050 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 26 | 0 | 0 | 37.1 |
| samtools | h020s | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |

### SNP distances

| Run | Method | Samples | Best convention | Exact pairs | Max error | Wrong pairs |
|---|---|---|---|---|---|---|
| flye | ska | 10/10 | events | 36/45 | 4 | Reference-h050: 4 (truth 0); h000-h050: 4 (truth 0); h002-h050: 4 (truth 0); h005-h050: 4 (truth 0); h010-h050: 4 (truth 0); h020-h050: 4 (truth 0) … |
| samtools | ska | 10/10 | events | 45/45 | 0 |  |

### Minor-haplotype sites

At the 30 SNP sites where the minor haplotype differs (inverted-repeat SNPs count twice), what the assembly has: the major allele (correct), the minor allele, N, another base, both alleles in different contigs or repeat copies (split) or nothing (missing); then the 3 indels (mdel1, mins1, mdel2): M major, m minor, N, o other, s split, - missing. Distance: SNPs to the reference reported by each method (truth).

| Run | Sample | Minor fraction | Major | Minor | N | Other | Split | Missing | Indels | N bases | Distance to reference (truth) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| flye | h000 | 0.00 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 0 (0) |
| flye | h002 | 0.02 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 0 (0) |
| flye | h005 | 0.05 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 0 (0) |
| flye | h010 | 0.10 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 0 (0) |
| flye | h020 | 0.20 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 0 (0) |
| flye | h020s | 0.20 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 12 (12) |
| flye | h030 | 0.30 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 0 (0) |
| flye | h040 | 0.40 | 24 (80%) | 0 (0%) | 0 (0%) | 0 | 6 | 0 | M/M/M | 0 | ska 0 (0) |
| flye | h050 | 0.50 | 18 (60%) | 4 (13%) | 0 (0%) | 0 | 8 | 0 | M/o/M | 0 | ska 4 (0) |
| samtools | h000 | 0.00 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 0 (0) |
| samtools | h002 | 0.02 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 0 (0) |
| samtools | h005 | 0.05 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 0 (0) |
| samtools | h010 | 0.10 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 0 (0) |
| samtools | h020 | 0.20 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 0 (0) |
| samtools | h020s | 0.20 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | ska 12 (12) |
| samtools | h030 | 0.30 | 25 (83%) | 0 (0%) | 5 (17%) | 0 | 0 | 0 | M/M/M | 5 | ska 0 (0) |
| samtools | h040 | 0.40 | 7 (23%) | 0 (0%) | 19 (63%) | 0 | 4 | 0 | M/M/M | 21 | ska 0 (0) |
| samtools | h050 | 0.50 | 0 (0%) | 0 (0%) | 22 (73%) | 0 | 8 | 0 | M/M/M | 25 | ska 0 (0) |

| Run | Wall time of the full run (s) |
|---|---|
| flye | 198.3 |
| samtools | 12.78 |
