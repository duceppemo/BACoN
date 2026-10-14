## Scenario `numt`

### Assemblies

| Run | Sample | Contigs | Circular (truth) | Length / truth | Covered | Duplicated bp | Mismatches | Ins bp | Del bp | Errors /100 kb |
|---|---|---|---|---|---|---|---|---|---|---|
| flye | n01 | 1 | 1 (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | n02 | 1 | 1 (1) | 120,001 / 120,000 | 1.0000 | 1 | 0 | 0 | 0 | 0.0 |
| flye | n03 | 3 | 1 (1) | 185,011 / 120,000 | 1.0000 | 11,630 | 5 | 2 | 0 | 5.3 |
| flye | n04 | 6 | 0 (1) | 396,615 / 120,000 | 1.0000 | 40,065 | 336 | 45 | 50 | 269.3 |
| flye | n05 | 1 | 1 (1) | 119,830 / 120,000 | 1.0000 | 0 | 1 | 0 | 0 | 0.8 |
| flye | n06 | 6 | 0 (1) | 376,768 / 120,000 | 1.0000 | 39,716 | 35 | 3 | 0 | 23.8 |
| samtools | n01 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | n02 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | n03 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 1 | 0 | 0 | 0.8 |
| samtools | n04 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 729 | 0 | 0 | 607.5 |
| samtools | n05 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | n06 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 32 | 0 | 0 | 26.7 |

### SNP distances

| Run | Method | Samples | Best convention | Exact pairs | Max error | Wrong pairs |
|---|---|---|---|---|---|---|
| flye | ska | 7/7 | events | 16/21 | 3 | Reference-n06: 4 (truth 7); n01-n06: 4 (truth 7); n02-n06: 7 (truth 10); n03-n06: 9 (truth 12); n05-n06: 7 (truth 10) |
| samtools | ska | 7/7 | events | 7/21 | 4 | Reference-n03: 4 (truth 5); Reference-n04: 2 (truth 5); Reference-n06: 4 (truth 7); n01-n03: 4 (truth 5); n01-n04: 2 (truth 5); n01-n06: 4 (truth 7) … |

### Reads by origin

Reads simulated / baited (1_extracted) / kept by the filter (2_filtered). NUPT reads: nuclear reads covering at least 100 bp of an insert; nuclear: the other nuclear reads. Share: of the filtered reads (bases).

| Run | Sample | Nuclear depth (genome) | Organelle reads | NUPT reads | Nuclear reads | NUPT share of filtered reads (bases) |
|---|---|---|---|---|---|---|
| flye | n01 | 0.5x (A) | 808/808/716 | 5/5/5 | 16/0/0 | 0.7% (0.8%) |
| flye | n02 | 2x (A) | 855/855/749 | 26/26/25 | 52/0/0 | 3.2% (4.7%) |
| flye | n03 | 10x (A) | 833/833/724 | 139/138/127 | 312/1/1 | 14.9% (17.9%) |
| flye | n04 | 25x (A) | 829/829/724 | 344/344/307 | 807/4/4 | 29.7% (33.6%) |
| flye | n05 | 10x (B) | 884/884/768 | 141/141/129 | 313/1/1 | 14.4% (16.8%) |
| flye | n06 | 25x (B) | 815/815/699 | 372/372/337 | 768/4/4 | 32.4% (35.1%) |
| samtools | n01 | 0.5x (A) | 808/808/716 | 5/5/5 | 16/0/0 | 0.7% (0.8%) |
| samtools | n02 | 2x (A) | 855/855/749 | 26/26/25 | 52/0/0 | 3.2% (4.7%) |
| samtools | n03 | 10x (A) | 833/833/724 | 139/138/127 | 312/1/1 | 14.9% (17.9%) |
| samtools | n04 | 25x (A) | 829/829/724 | 344/344/307 | 807/4/4 | 29.7% (33.6%) |
| samtools | n05 | 10x (B) | 884/884/768 | 141/141/129 | 313/1/1 | 14.4% (16.8%) |
| samtools | n06 | 25x (B) | 815/815/699 | 372/372/337 | 768/4/4 | 32.4% (35.1%) |

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
| flye | n03 | 0 (45,000) | 7 (86,630) | 1703 / 0 / 0 / 0 / 0 / 0 | 0 | 2 / 0 / 0 / 0 | ska 5 (5) |
| flye | n04 | 424 (74,424) | 7 (85,641) | 1370 / 0 / 0 / 0 / 332 / 0 | 332 | 3 / 0 / 0 / 0 | ska 5 (5) |
| flye | n05 | 0 (45,000) | 1 (74,830) | 46 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | ska 3 (3) |
| flye | n06 | 34 (74,767) | 4 (84,949) | 15 / 0 / 0 / 0 / 31 / 0 | 31 | 0 / 0 / 0 / 3 | ska 4 (7) |
| samtools | n01 | 0 (45,000) | 0 (75,000) | 1703 / 0 / 0 / 0 / 0 / 0 | 0 | 0 / 0 / 0 / 0 | ska 0 (0) |
| samtools | n02 | 0 (45,000) | 0 (75,000) | 1703 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | ska 3 (3) |
| samtools | n03 | 1 (45,000) | 0 (75,000) | 1702 / 0 / 1 / 0 / 0 / 0 | 0 | 2 / 0 / 0 / 0 | ska 4 (5) |
| samtools | n04 | 729 (45,000) | 0 (75,000) | 981 / 0 / 721 / 0 / 0 / 0 | 0 | 3 / 0 / 0 / 0 | ska 2 (5) |
| samtools | n05 | 0 (45,000) | 0 (75,000) | 46 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | ska 3 (3) |
| samtools | n06 | 32 (45,000) | 0 (75,000) | 17 / 0 / 29 / 0 / 0 / 0 | 0 | 0 / 0 / 3 / 0 | ska 4 (7) |

| Run | Wall time of the full run (s) |
|---|---|
| flye | 109.01 |
| samtools | 8.0 |
