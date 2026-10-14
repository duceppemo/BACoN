## Scenario `numt`

### Assemblies

| Run | Sample | Contigs | Circular (truth) | Length / truth | Covered | Duplicated bp | Mismatches | Ins bp | Del bp | Errors /100 kb |
|---|---|---|---|---|---|---|---|---|---|---|
| flye | n01 | 1 | 1 (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | n02 | 1 | 1 (1) | 120,001 / 120,000 | 1.0000 | 1 | 0 | 0 | 0 | 0.0 |
| flye | n03 | 2 | 1 (1) | 149,233 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | n04 | 7 | 2 (1) | 339,412 / 120,000 | 1.0000 | 42,931 | 190 | 35 | 24 | 152.8 |
| flye | n05 | 1 | 1 (1) | 119,830 / 120,000 | 1.0000 | 0 | 1 | 0 | 0 | 0.8 |
| flye | n06 | 7 | 0 (1) | 183,726 / 120,000 | 0.9999 | 12,879 | 19 | 8 | 0 | 20.3 |
| samtools | n01 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | n02 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | n03 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 5 | 0 | 0 | 4.2 |
| samtools | n04 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 756 | 0 | 0 | 630.0 |
| samtools | n05 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | n06 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 30 | 0 | 0 | 25.0 |

### SNP distances

| Run | Method | Samples | Best convention | Exact pairs | Max error | Wrong pairs |
|---|---|---|---|---|---|---|
| flye | ska | 7/7 | events | 21/21 | 0 |  |
| samtools | ska | 7/7 | events | 7/21 | 4 | Reference-n03: 4 (truth 5); Reference-n04: 2 (truth 5); Reference-n06: 4 (truth 7); n01-n03: 4 (truth 5); n01-n04: 2 (truth 5); n01-n06: 4 (truth 7) … |

### Reads by origin

Reads simulated / baited (1_extracted) / kept by the filter (2_filtered). NUPT reads: nuclear reads covering at least 100 bp of an insert; nuclear: the other nuclear reads. Share: of the filtered reads (bases).

| Run | Sample | Nuclear depth (genome) | Organelle reads | NUPT reads | Nuclear reads | NUPT share of filtered reads (bases) |
|---|---|---|---|---|---|---|
| flye | n01 | 0.5x (A) | 808/808/707 | 5/5/0 | 16/0/0 | 0.0% (0.0%) |
| flye | n02 | 2x (A) | 855/855/825 | 26/26/3 | 52/0/0 | 0.4% (0.4%) |
| flye | n03 | 10x (A) | 833/833/833 | 139/138/107 | 312/1/0 | 11.4% (13.2%) |
| flye | n04 | 25x (A) | 829/829/829 | 344/344/302 | 807/4/0 | 26.7% (30.2%) |
| flye | n05 | 10x (B) | 884/884/884 | 141/141/106 | 313/1/0 | 10.7% (12.2%) |
| flye | n06 | 25x (B) | 815/815/815 | 372/372/336 | 768/4/0 | 29.2% (31.7%) |
| samtools | n01 | 0.5x (A) | 808/808/707 | 5/5/0 | 16/0/0 | 0.0% (0.0%) |
| samtools | n02 | 2x (A) | 855/855/825 | 26/26/3 | 52/0/0 | 0.4% (0.4%) |
| samtools | n03 | 10x (A) | 833/833/833 | 139/138/107 | 312/1/0 | 11.4% (13.2%) |
| samtools | n04 | 25x (A) | 829/829/829 | 344/344/302 | 807/4/0 | 26.7% (30.2%) |
| samtools | n05 | 10x (B) | 884/884/884 | 141/141/106 | 313/1/0 | 10.7% (12.2%) |
| samtools | n06 | 25x (B) | 815/815/815 | 372/372/336 | 768/4/0 | 29.2% (31.7%) |

### NUPT reads baited, by insert

Baited / simulated reads of each insert (length, identity to the plastome in genome A; in genome B every insert is 99.9% identical).

| Run | Sample | #1 0.5 kb 99.9% | #2 0.5 kb 90% | #3 2 kb 99% | #4 2 kb 97% | #5 5 kb 97% | #6 5 kb 99.9% | #7 15 kb 90% | #8 15 kb 99% |
|---|---|---|---|---|---|---|---|---|---|
| flye | n01 | 0/0 | 0/0 | 1/1 | 0/0 | 1/1 | 1/1 | 2/2 | 0/0 |
| flye | n02 | 3/3 | 1/1 | 2/2 | 3/3 | 0/0 | 2/2 | 4/4 | 11/11 |
| flye | n03 | 8/8 | 14/15 | 5/5 | 9/9 | 19/19 | 17/17 | 37/37 | 29/29 |
| flye | n04 | 24/24 | 28/28 | 23/23 | 21/21 | 42/42 | 44/44 | 72/72 | 90/90 |
| flye | n05 | 10/10 | 14/14 | 12/12 | 13/13 | 12/12 | 17/17 | 39/39 | 24/24 |
| flye | n06 | 20/20 | 26/26 | 31/31 | 27/27 | 50/50 | 47/47 | 73/73 | 98/98 |
| samtools | n01 | 0/0 | 0/0 | 1/1 | 0/0 | 1/1 | 1/1 | 2/2 | 0/0 |
| samtools | n02 | 3/3 | 1/1 | 2/2 | 3/3 | 0/0 | 2/2 | 4/4 | 11/11 |
| samtools | n03 | 8/8 | 14/15 | 5/5 | 9/9 | 19/19 | 17/17 | 37/37 | 29/29 |
| samtools | n04 | 24/24 | 28/28 | 23/23 | 21/21 | 42/42 | 44/44 | 72/72 | 90/90 |
| samtools | n05 | 10/10 | 14/14 | 12/12 | 13/13 | 12/12 | 17/17 | 39/39 | 24/24 |
| samtools | n06 | 20/20 | 26/26 | 31/31 | 27/27 | 50/50 | 47/47 | 73/73 | 98/98 |

### Consensus

Consensus errors (mismatched, inserted and deleted bp) inside the plastome segments copied into the nucleus and elsewhere, per aligned length; at the substitutions that distinguish the NUPTs from the plastome, what the assembly has (plastome base / NUPT base / N / other / split: different bases in different contigs / missing) and the number of sites where the NUPT base is in any contig; the sample's own SNPs inside the copied segments (where the NUPT reads carry the reference base): called / reference base / N / split or other; distance: SNPs to the reference reported by each method (truth).

| Run | Sample | Errors in copied segments (aligned bp) | Errors elsewhere (aligned bp) | NUPT sites: plastome / NUPT / N / other / split / missing | NUPT base in a contig | Own SNPs in copied segments | Distance to reference (truth) |
|---|---|---|---|---|---|---|---|
| flye | n01 | 0 (45,000) | 0 (75,000) | 1703 / 0 / 0 / 0 / 0 / 0 | 0 | 0 / 0 / 0 / 0 | ska 0 (0) |
| flye | n02 | 0 (45,000) | 0 (75,001) | 1703 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | ska 3 (3) |
| flye | n03 | 0 (45,000) | 0 (75,000) | 1703 / 0 / 0 / 0 / 0 / 0 | 0 | 2 / 0 / 0 / 0 | ska 5 (5) |
| flye | n04 | 245 (72,435) | 4 (90,496) | 1551 / 0 / 0 / 0 / 151 / 0 | 151 | 3 / 0 / 0 / 0 | ska 5 (5) |
| flye | n05 | 0 (45,000) | 1 (74,830) | 46 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | ska 3 (3) |
| flye | n06 | 2 (46,992) | 25 (85,879) | 44 / 0 / 0 / 0 / 2 / 0 | 2 | 3 / 0 / 0 / 0 | ska 7 (7) |
| samtools | n01 | 0 (45,000) | 0 (75,000) | 1703 / 0 / 0 / 0 / 0 / 0 | 0 | 0 / 0 / 0 / 0 | ska 0 (0) |
| samtools | n02 | 0 (45,000) | 0 (75,000) | 1703 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | ska 3 (3) |
| samtools | n03 | 5 (45,000) | 0 (75,000) | 1698 / 0 / 5 / 0 / 0 / 0 | 0 | 2 / 0 / 0 / 0 | ska 4 (5) |
| samtools | n04 | 756 (45,000) | 0 (75,000) | 953 / 0 / 749 / 0 / 0 / 0 | 0 | 3 / 0 / 0 / 0 | ska 2 (5) |
| samtools | n05 | 0 (45,000) | 0 (75,000) | 46 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | ska 3 (3) |
| samtools | n06 | 30 (45,000) | 0 (75,000) | 19 / 0 / 27 / 0 / 0 / 0 | 0 | 0 / 0 / 3 / 0 | ska 4 (7) |

| Run | Wall time of the full run (s) |
|---|---|
| flye | 115.41 |
| samtools | 9.45 |
