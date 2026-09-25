## Scenario `rdna`

### Assemblies

| Run | Sample | Contigs | Circular (truth) | Length / truth | Covered | Duplicated bp | Mismatches | Ins bp | Del bp | Errors /100 kb |
|---|---|---|---|---|---|---|---|---|---|---|
| flye | t01 | 1 | 1 (0) | 6,000 / 22,000 | 0.2727 | 0 | 0 | 0 | 0 | 0.0 |
| flye | t02 | 1 | 1 (0) | 6,000 / 22,000 | 0.2727 | 0 | 0 | 0 | 0 | 0.0 |
| flye | t03 | 1 | 1 (0) | 6,000 / 22,000 | 0.2727 | 0 | 0 | 0 | 0 | 0.0 |
| flye | t04 | 1 | 1 (0) | 6,000 / 22,000 | 0.2727 | 0 | 0 | 0 | 0 | 0.0 |
| myloasm | t01 | 1 | 0 (0) | 22,000 / 22,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| myloasm | t02 | 1 | 0 (0) | 18,981 / 22,000 | 0.9991 | 0 | 0 | 0 | 3000 | 15805.3 |
| myloasm | t03 | 1 | 0 (0) | 21,972 / 22,000 | 0.9987 | 0 | 1 | 0 | 0 | 4.6 |
| myloasm | t04 | 1 | 0 (0) | 22,000 / 22,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | t01 | 1 | NA (0) | 22,000 / 22,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | t02 | 1 | NA (0) | 22,000 / 22,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | t03 | 1 | NA (0) | 22,000 / 22,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | t04 | 1 | NA (0) | 22,000 / 22,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |

### SNP distances

| Run | Method | Samples | Best convention | Exact pairs | Max error | Wrong pairs |
|---|---|---|---|---|---|---|
| flye | ska | 5/5 | no_repeat | 1/10 | 8 | Reference-t02: 1 (truth 3); Reference-t03: 1 (truth 5); Reference-t04: 1 (truth 5); t01-t02: 1 (truth 3); t01-t03: 1 (truth 5); t01-t04: 1 (truth 5) … |
| flye | ska_0.5 | 5/5 | no_repeat | 1/10 | 8 | Reference-t02: 1 (truth 3); Reference-t03: 1 (truth 5); Reference-t04: 1 (truth 5); t01-t02: 1 (truth 3); t01-t03: 1 (truth 5); t01-t04: 1 (truth 5) … |
| myloasm | ska | 5/5 | events | 6/10 | 1 | Reference-t03: 7 (truth 6); t01-t03: 7 (truth 6); t02-t03: 3 (truth 2); t03-t04: 13 (truth 12) |
| samtools | ska | 5/5 | events | 10/10 | 0 |  |

| Run | Wall time of the full run (s) |
|---|---|
| flye | 29.18 |
| myloasm | 3.33 |
| samtools | 1.55 |
