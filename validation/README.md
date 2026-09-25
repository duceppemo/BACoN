# Validation

Simulated datasets with a known truth, the scripts that run BACoN on them and score the results, and the
dated records of every validation.

| File | Purpose |
|---|---|
| `simulate.py` | Generates the scenarios (reference, reads, true genomes, `truth.json`). Deterministic: the same seed gives byte-identical files. Standard library only |
| `run_validation.sh` | Runs BACoN on every scenario with every assembler and SNP method, then `score.py` |
| `score.py` | Scores BACoN output folders against the truth: assembly completeness and consensus errors (minimap2), SNP distances |
| `compare_runs.py` | For real data without a truth: agreement between runs (assemblers, SNP methods) |
| `results/<date>_<version>/` | Records: summary, scores, program versions |

## Scenarios

| Scenario | Genome | Samples | Tests |
|---|---|---|---|
| `plastid` | 70 kb circular, chloroplast-like: large single copy, inverted repeat (10 kb), small single copy, inverted repeat | 9: two identical to the reference, two 1–2 SNPs away, a clade of three, one with 40 SNPs and a 1.5 kb deletion, one at 10x | identical and very close samples, SNPs inside the inverted repeat, a 5 bp insertion, a deletion, low depth |
| `linear` | three linear molecules (15, 9 and 6 kb) | 4 | assemblers circularizing or dropping short linear molecules |
| `rdna` | linear 22 kb: 2 kb flank, a 3 kb unit repeated 6 times, 2 kb flank | 4, with SNPs in the flanks and in every copy of the unit | tandem arrays collapsing |

Reads mimic R10.4.1 SUP reads: about 0.4% substitutions, 0.2% random insertions and deletions, and homopolymer
length errors that grow with the run length; log-normal lengths (N50 about 8 kb); 15% of the reads come from
an unrelated sequence, which baiting must remove. Real Nanopore errors are more systematic than these; the
real-data comparisons in the records complement the simulations.

SNP distances are compared with the truth under three conventions: `events` (each mutation once; a SNP in an
inverted repeat or in every copy of a tandem unit is one event), `sites` (each changed position) and
`no_repeat` (SNPs outside repeats only). The convention that fits a method best says how it treats repeats.

## Running

```bash
conda activate BACoN
bash validation/run_validation.sh /tmp/bacon_validation 16
```

It takes a few minutes. To add a record, copy `versions.txt` and `SCORES_*` into a new
`results/<date>_<version>/` folder with a `SUMMARY.md`. Records are never rewritten: an error found later is
noted as an erratum in a new record.
