# Validation of BACoN 0.3.8 (2026-10-13)

A rerun of the validation of 0.3.0 ([`../2026-09-24_v0.3.0/`](../2026-09-24_v0.3.0/SUMMARY.md)) with the released
0.3.8 and the same program versions, after the fixes of 0.3.1 to 0.3.8. The screening of candidate programs
was not repeated: it decided the programs, and none was added since.

All runs (simulated and potato) used one conda environment, with BACoN 0.3.8 installed from the v0.3.8 tag, on
one machine: 2 Intel Xeon Gold 6226R (64 threads in all), 503 GB of RAM, Linux.

## Simulated data

```bash
bash validation/run_validation.sh <folder> 16      # BACoN 0.3.8 installed from the v0.3.8 tag
```

Program versions: [`versions.txt`](versions.txt) (minimap2 2.31, Filtlong 0.3.1, Flye 2.9.6, myloasm 0.7.0,
samtools 1.24, SKA2 0.5.1, Parsnp 2.1.5, HarvestTools 1.3, FastTree 2.2.0: those of 0.3.0). Scores:
[`SCORES_plastid.md`](SCORES_plastid.md), [`SCORES_linear.md`](SCORES_linear.md), [`SCORES_rdna.md`](SCORES_rdna.md)
(and `.json`); wall times: [`timings.txt`](timings.txt).

**Every score is identical to 0.3.0**: assemblies (samples, contigs, circular contigs, consensus errors,
completeness) and SNP distances, for every scenario, assembler and SNP method; the score files differ only in
their wall times. The tables of the 0.3.0 summary therefore hold for 0.3.8:

| Scenario | samtools + SKA2 (default) | Flye + SKA2 | myloasm + SKA2 | Flye + SKA2 pan (0.5) | Flye + Parsnp `-c` | Flye + Parsnp default |
|---|---|---|---|---|---|---|
| plastid | **45/45** | 45/45 | 45/45 | 45/45 | 36/45 | 22/36, drops the 10x sample |
| linear | **10/10** | 10/10 | 6/10 | 10/10 | 10/10 | 10/10 |
| rdna | **10/10** | 1/10 | 6/10 | 1/10 | fails (assemblies 3.7 times shorter than the reference) | — |

Pairs of genomes with exactly the true SNP distance. The failure of Parsnp on the rdna scenario is the
expected one: Flye collapses each tandem array into one circular unit, and Parsnp rejects the assemblies.

Wall time, plastid scenario (9 samples, `-t 16 -p 4`): samtools 4.4 s, myloasm 10.0 s, Flye 58.8 s
(0.3.0: 6.9, 11.8 and 60.3 s). Peak memory: about 240 MB (samtools, Flye) and 2.0 GB (myloasm).

The assemblies themselves were compared with those of the 0.3.0 record: the 17 templated and the 17 myloasm
assemblies are byte-identical (myloasm, found not fully deterministic in 0.3.0, was not this time); 16 of the 17
Flye assemblies are, and one (linear, `l03`) differs in one contig, with the same score. **Flye is not fully
deterministic either.**

Changes since 0.3.0 that these scores do not exercise: the number of SNP sites (0.3.5 counts only columns with
two or more nucleotides; distances are unchanged), annotations and SNP effects (0.3.6), sample metadata, the
report, Parsnp 2.1.1's off-by-one positions (0.3.7; this record uses Parsnp 2.1.5), and IQ-TREE (the scores use
FastTree). Running the test suite the same day found that `--tree iqtree` fails in 0.3.8 when IQ-TREE keeps
fewer than four sequences (it keeps two of each set of identical sequences, then refuses a bootstrap); fixed
after 0.3.8 (CHANGELOG, "Unreleased").

## Real data (public)

The 28 potato cultivars of BioProject PRJNA807056 (MinION, R9.4.1-era reads, mean Q14–16, while the templated
consensus uses the R10.4 SUP error profile and the simulations mimic R10.4.1 reads; reference NC_008096.2), as in the
[tutorial](https://github.com/duceppemo/BACoN/wiki/Tutorial), `-t 48 -p 12`, compared with
`validation/compare_runs.py` ([`potato_compare_runs.txt`](potato_compare_runs.txt), times in
[`potato_timings.txt`](potato_timings.txt)):

| Comparison | Pairs with identical SNP distances (of 406) | Largest difference | 0.3.0 |
|---|---|---|---|
| samtools + SKA2 vs Flye + SKA2 | **406** | 0 | 406, 0 |
| samtools + SKA2 vs Flye + Parsnp `-c` | 165 | 11 | 121, 10 (vs Flye + SKA2) |
| samtools + SKA2 vs myloasm + SKA2 | 93 | 21 | 93, 21 (vs Flye + SKA2) |

The templated and the Flye assemblies, analysed independently, still agree exactly. Flye assembled 17 of the
28 plastomes as one circular contig (as in 0.3.0); myloasm gave a circular contig in 8 samples, a single
contig in 2 of them; every templated assembly is within 0.1% of the
reference length. SNP sites: 118 (templated and Flye, SKA2), 127 (Flye, Parsnp), 93 (myloasm). Wall time on
28 samples: templated 124 s, myloasm 151 s, Flye 725 s (peak memory 0.7, 2.1 and 0.7 GB).

Two of the 28 Flye assemblies (`15_22_4`, `Bankir`) differ from those of the 0.3.0 run, the other 26 are
byte-identical: Flye's run-to-run variation, as on the simulated data. The SKA2 distances of the two Flye runs
are identical, but Parsnp's alignment of the 28 assemblies changed (the same 127 sites; 78 columns differ), which
is why Flye + Parsnp agrees on 165 pairs here and on 121 in 0.3.0: a change in the input assemblies, not in
BACoN's use of Parsnp (Parsnp 2.1.5 in both). Parsnp's distances move with small changes of the assemblies;
SKA2's did not.
