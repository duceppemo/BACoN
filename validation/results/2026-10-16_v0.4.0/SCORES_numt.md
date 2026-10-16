## Scenario `numt`

### Assemblies

| Run | Sample | Contigs | Circular (truth) | Length / truth | Covered | Duplicated bp | Mismatches | Ins bp | Del bp | Errors /100 kb |
|---|---|---|---|---|---|---|---|---|---|---|
| flye | n01 | 1 | 1 (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | n02 | 1 | 1 (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | n03 | 1 | 1 (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | n04 | 6 | 1 (1) | 276,082 / 120,000 | 1.0000 | 4,783 | 7 | 3 | 0 | 8.0 |
| flye | n05 | 1 | 1 (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| flye | n06 | 9 | 1 (1) | 313,104 / 120,000 | 1.0000 | 20,733 | 16 | 2 | 0 | 12.8 |
| samtools | n01 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | n02 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | n03 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 2 | 0 | 0 | 1.7 |
| samtools | n04 | 1 | NA (1) | 120,003 / 120,000 | 1.0000 | 0 | 1230 | 3 | 0 | 1027.5 |
| samtools | n05 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | n06 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 22 | 0 | 0 | 18.3 |

### SNP distances

| Run | Method | Samples | Best convention | Exact pairs | Max error | Wrong pairs |
|---|---|---|---|---|---|---|
| flye | minimap2 | 7/7 | events | 10/21 | 3 | Reference-n04: 4 (truth 5); Reference-n06: 4 (truth 7); n01-n04: 4 (truth 5); n01-n06: 4 (truth 7); n02-n04: 7 (truth 8); n02-n06: 7 (truth 10) … |
| flye | ska | 7/7 | events | 16/21 | 2 | Reference-n06: 5 (truth 7); n01-n06: 5 (truth 7); n02-n06: 8 (truth 10); n03-n06: 10 (truth 12); n05-n06: 8 (truth 10) |
| samtools | minimap2 | 7/7 | events | 11/21 | 3 | Reference-n04: 2 (truth 5); Reference-n06: 4 (truth 7); n01-n04: 2 (truth 5); n01-n06: 4 (truth 7); n02-n04: 5 (truth 8); n02-n06: 7 (truth 10) … |
| samtools | ska | 7/7 | events | 11/21 | 3 | Reference-n04: 2 (truth 5); Reference-n06: 4 (truth 7); n01-n04: 2 (truth 5); n01-n06: 4 (truth 7); n02-n04: 5 (truth 8); n02-n06: 7 (truth 10) … |

### Reads by origin

Reads simulated / baited (1_extracted) / kept by the filter (2_filtered). NUPT reads: nuclear reads covering at least 100 bp of an insert; nuclear: the other nuclear reads. Share: of the filtered reads (bases).

| Run | Sample | Nuclear depth (genome) | Organelle reads | NUPT reads | Nuclear reads | NUPT share of filtered reads (bases) |
|---|---|---|---|---|---|---|
| flye | n01 | 0.5x (A) | 838/838/742 | 7/7/6 | 48/1/1 | 0.8% (1.4%) |
| flye | n02 | 2x (A) | 839/839/745 | 28/28/24 | 190/0/0 | 3.1% (3.7%) |
| flye | n03 | 10x (A) | 819/819/724 | 146/144/133 | 935/2/2 | 15.5% (17.0%) |
| flye | n04 | 25x (A) | 809/809/694 | 341/339/303 | 2361/0/0 | 30.4% (33.6%) |
| flye | n05 | 10x (B) | 829/829/738 | 143/143/132 | 951/1/1 | 15.2% (17.7%) |
| flye | n06 | 25x (B) | 845/845/731 | 343/343/308 | 2286/0/0 | 29.6% (34.8%) |
| samtools | n01 | 0.5x (A) | 838/838/742 | 7/7/6 | 48/1/1 | 0.8% (1.4%) |
| samtools | n02 | 2x (A) | 839/839/745 | 28/28/24 | 190/0/0 | 3.1% (3.7%) |
| samtools | n03 | 10x (A) | 819/819/724 | 146/144/133 | 935/2/2 | 15.5% (17.0%) |
| samtools | n04 | 25x (A) | 809/809/694 | 341/339/303 | 2361/0/0 | 30.4% (33.6%) |
| samtools | n05 | 10x (B) | 829/829/738 | 143/143/132 | 951/1/1 | 15.2% (17.7%) |
| samtools | n06 | 25x (B) | 845/845/731 | 343/343/308 | 2286/0/0 | 29.6% (34.8%) |

### NUPT reads baited, by insert

Baited / simulated reads of each insert (length, identity to the plastome in genome A; in genome B every insert is 99.9% identical).

| Run | Sample | #1 0.5 kb 99.9% | #2 0.5 kb 90% | #3 2 kb 99% | #4 2 kb 97% | #5 5 kb 97% | #6 5 kb 99.9% | #7 15 kb 90% | #8 15 kb 99% |
|---|---|---|---|---|---|---|---|---|---|
| flye | n01 | 0/0 | 0/0 | 0/0 | 1/1 | 2/2 | 1/1 | 1/1 | 2/2 |
| flye | n02 | 1/1 | 0/0 | 2/2 | 2/2 | 4/4 | 5/5 | 6/6 | 8/8 |
| flye | n03 | 10/10 | 12/12 | 18/18 | 14/14 | 13/13 | 12/12 | 28/30 | 37/37 |
| flye | n04 | 21/21 | 21/22 | 19/19 | 25/25 | 29/29 | 47/47 | 88/89 | 89/89 |
| flye | n05 | 10/10 | 13/13 | 7/7 | 10/10 | 27/27 | 15/15 | 33/33 | 28/28 |
| flye | n06 | 28/28 | 30/30 | 23/23 | 30/30 | 47/47 | 45/45 | 78/78 | 62/62 |
| samtools | n01 | 0/0 | 0/0 | 0/0 | 1/1 | 2/2 | 1/1 | 1/1 | 2/2 |
| samtools | n02 | 1/1 | 0/0 | 2/2 | 2/2 | 4/4 | 5/5 | 6/6 | 8/8 |
| samtools | n03 | 10/10 | 12/12 | 18/18 | 14/14 | 13/13 | 12/12 | 28/30 | 37/37 |
| samtools | n04 | 21/21 | 21/22 | 19/19 | 25/25 | 29/29 | 47/47 | 88/89 | 89/89 |
| samtools | n05 | 10/10 | 13/13 | 7/7 | 10/10 | 27/27 | 15/15 | 33/33 | 28/28 |
| samtools | n06 | 28/28 | 30/30 | 23/23 | 30/30 | 47/47 | 45/45 | 78/78 | 62/62 |

### Consensus

Consensus errors (mismatched, inserted and deleted bp) inside the plastome segments copied into the nucleus and elsewhere, per aligned length; at the substitutions that distinguish the NUPTs from the plastome, what the assembly has (plastome base / NUPT base / N / other / split: different bases in different contigs / missing) and the number of sites where the NUPT base is in any contig; the sample's own SNPs inside the copied segments (where the NUPT reads carry the reference base): called / reference base / N / split or other; distance: SNPs to the reference reported by each method (truth).

| Run | Sample | Errors in copied segments (aligned bp) | Errors elsewhere (aligned bp) | NUPT sites: plastome / NUPT / N / other / split / missing | NUPT base in a contig | Own SNPs in copied segments | Distance to reference (truth) |
|---|---|---|---|---|---|---|---|
| flye | n01 | 0 (45,000) | 0 (75,000) | 1666 / 0 / 0 / 0 / 0 / 0 | 0 | 0 / 0 / 0 / 0 | minimap2 0; ska 0 (0) |
| flye | n02 | 0 (45,000) | 0 (75,000) | 1666 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | minimap2 3; ska 3 (3) |
| flye | n03 | 0 (45,000) | 0 (75,000) | 1666 / 0 / 0 / 0 / 0 / 0 | 0 | 2 / 0 / 0 / 0 | minimap2 5; ska 5 (5) |
| flye | n04 | 5 (50,489) | 5 (74,294) | 1661 / 0 / 0 / 0 / 5 / 0 | 5 | 3 / 0 / 0 / 0 | minimap2 4; ska 5 (5) |
| flye | n05 | 0 (45,000) | 0 (75,000) | 46 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | minimap2 3; ska 3 (3) |
| flye | n06 | 12 (54,500) | 6 (86,233) | 36 / 0 / 0 / 0 / 10 / 0 | 10 | 1 / 0 / 0 / 2 | minimap2 4; ska 5 (7) |
| samtools | n01 | 0 (45,000) | 0 (75,000) | 1666 / 0 / 0 / 0 / 0 / 0 | 0 | 0 / 0 / 0 / 0 | minimap2 0; ska 0 (0) |
| samtools | n02 | 0 (45,000) | 0 (75,000) | 1666 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | minimap2 3; ska 3 (3) |
| samtools | n03 | 2 (45,000) | 0 (75,000) | 1664 / 0 / 2 / 0 / 0 / 0 | 0 | 2 / 0 / 0 / 0 | minimap2 5; ska 5 (5) |
| samtools | n04 | 1233 (45,000) | 0 (75,000) | 447 / 0 / 1219 / 0 / 0 / 0 | 0 | 2 / 0 / 1 / 0 | minimap2 2; ska 2 (5) |
| samtools | n05 | 0 (45,000) | 0 (75,000) | 46 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | minimap2 3; ska 3 (3) |
| samtools | n06 | 22 (45,000) | 0 (75,000) | 26 / 0 / 20 / 0 / 0 / 0 | 0 | 1 / 0 / 2 / 0 | minimap2 4; ska 4 (7) |

| Run | Wall time of the full run (s) |
|---|---|
| flye | 130.45 |
| samtools | 8.8 |
