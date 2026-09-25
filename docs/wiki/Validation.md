# Validation

BACoN 0.3 was rebuilt around one rule: a program is used only if it is maintained and does at least as well
as the alternatives on data where the answer is known. The records, with every score, are in
[`validation/results/`](https://github.com/duceppemo/BACoN/tree/main/validation/results); this page
summarizes the [2026-09-24 record](https://github.com/duceppemo/BACoN/blob/main/validation/results/2026-09-24_v0.3.0/SUMMARY.md).

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
| Plastid sample at 10x | 1 contig, 1 error | 3 contigs, 68 errors | 2 contigs, 0 errors |
| Linear molecules | correct | correct, none circularized | correct |
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
| PhaME | its dependency check fails with any samtools from 1.10 |

## Real data

The 28 potato cultivars of the [Tutorial](Tutorial) (public, R9.4.1-era reads): the templated assembly and the
Flye assembly, analysed independently, give **identical SNP distances for all 406 pairs** of cultivars. They
recover the known cytoplasm types, including the 241 bp deletion that marks the T-type (resolved exactly by
Flye).

## Reproducing

```bash
conda activate BACoN
bash validation/run_validation.sh /tmp/bacon_validation 16
```

The simulated reads are generated from a fixed seed (byte-identical files). All results reproduce exactly,
except that myloasm is not fully deterministic: its assembly of the 10x sample differed by one base between
two runs.
