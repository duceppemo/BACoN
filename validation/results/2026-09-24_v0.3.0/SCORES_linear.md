## Scenario `linear`

### Assemblies

| Run | Sample | Contigs | Circular (truth) | Length / truth | Covered | Duplicated bp | Mismatches | Ins bp | Del bp | Errors /100 kb |
|---|---|---|---|---|---|---|---|---|---|---|
| flye | l01 | 3 | 0 (0) | 29,997 / 30,000 | 0.9999 | 0 | 0 | 0 | 0 | 0.0 |
| flye | l02 | 3 | 0 (0) | 29,223 / 30,000 | 0.9741 | 0 | 0 | 0 | 0 | 0.0 |
| flye | l03 | 3 | 0 (0) | 29,997 / 30,000 | 0.9999 | 0 | 0 | 0 | 0 | 0.0 |
| flye | l04 | 3 | 0 (0) | 28,780 / 30,000 | 0.9593 | 0 | 0 | 0 | 0 | 0.0 |
| myloasm | l01 | 3 | 0 (0) | 29,620 / 30,000 | 0.9873 | 0 | 2 | 0 | 0 | 6.8 |
| myloasm | l02 | 3 | 0 (0) | 29,892 / 30,000 | 0.9962 | 0 | 0 | 0 | 0 | 0.0 |
| myloasm | l03 | 3 | 0 (0) | 29,947 / 30,000 | 0.9982 | 0 | 0 | 0 | 0 | 0.0 |
| myloasm | l04 | 3 | 0 (0) | 29,763 / 30,000 | 0.9914 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | l01 | 3 | NA (0) | 30,000 / 30,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | l02 | 3 | NA (0) | 30,000 / 30,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | l03 | 3 | NA (0) | 30,000 / 30,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | l04 | 3 | NA (0) | 30,000 / 30,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |

### SNP distances

| Run | Method | Samples | Best convention | Exact pairs | Max error | Wrong pairs |
|---|---|---|---|---|---|---|
| flye | parsnp | 5/5 | events | 10/10 | 0 |  |
| flye | parsnp_default | 5/5 | events | 10/10 | 0 |  |
| flye | ska | 5/5 | events | 10/10 | 0 |  |
| flye | ska_0.5 | 5/5 | events | 10/10 | 0 |  |
| myloasm | ska | 5/5 | events | 6/10 | 2 | Reference-l01: 2 (truth 0); l01-l02: 7 (truth 5); l01-l03: 10 (truth 8); l01-l04: 14 (truth 12) |
| samtools | ska | 5/5 | events | 10/10 | 0 |  |

| Run | Wall time of the full run (s) |
|---|---|
| flye | 11.68 |
| myloasm | 3.16 |
| samtools | 1.5 |
