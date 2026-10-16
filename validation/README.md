# Validation

Simulated datasets with a known truth, the scripts that run BACoN on them and score the results, and the
dated records of every validation.

| File | Purpose |
|---|---|
| `simulate.py` | Generates the scenarios (reference, reads, true genomes, `truth.json`). Deterministic: the same seed gives byte-identical files. Standard library only |
| `run_validation.sh` | Runs BACoN on every scenario with every assembler and SNP method, then `score.py` (`new`: the heteroplasmy and NUPT scenarios, templated and Flye, compared with SKA2 and minimap2) |
| `score.py` | Scores BACoN output folders against the truth: assembly completeness and consensus errors (minimap2), SNP distances; for the new scenarios, the alleles called at minor-haplotype or NUPT sites and the origin of the baited reads |
| `compare_runs.py` | For real data without a truth: agreement between runs (assemblers, SNP methods) |
| `results/<date>_<version>/` | Records: summary, scores, program versions |

## Scenarios

| Scenario | Genome | Samples | Tests |
|---|---|---|---|
| `plastid` | 70 kb circular, chloroplast-like: large single copy, inverted repeat (10 kb), small single copy, inverted repeat | 9: two identical to the reference, two 1–2 SNPs away, a clade of three, one with 40 SNPs and a 1.5 kb deletion, one at 10x | identical and very close samples, SNPs inside the inverted repeat, a 5 bp insertion, a deletion, low depth |
| `linear` | three linear molecules (15, 9 and 6 kb) | 4 | assemblers circularizing or dropping short linear molecules |
| `rdna` | linear 22 kb: 2 kb flank, a 3 kb unit repeated 6 times, 2 kb flank | 4, with SNPs in the flanks and in every copy of the unit | tandem arrays collapsing |
| `heteroplasmy` | 70 kb circular, as `plastid` | 9 mixtures of a major and a minor haplotype (25 SNPs: 15 LSC, 5 SSC, 5 in the inverted repeat; and 3 indels of 1–5 bp), minor fraction 0, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40 and 0.50 (a tie); one more at 0.20 whose haplotypes both carry 12 SNPs of their own; 80x | at which minor fraction the consensus stops reporting the major allele (N, or the minor allele), and whether the SNP distances move |
| `numt` | 120 kb circular, chloroplast-like (LSC 80 kb, IR 12 kb, SSC 16 kb), 50x | 6 (0–7 SNPs), whose off-target reads are nuclear: four with reads of a 765 kb nuclear sequence holding eight plastome segments (NUPTs: 0.5, 2, 5 and 15 kb at 99.9, 99, 97 and 90% identity) between 80 kb spacers, at 0.5, 2, 10 and 25x; two (10x, 25x) whose NUPTs are all 99.9% identical, with substitutions recorded in the truth | how many NUPT reads baiting keeps, consensus errors inside the copied segments, whether NUPT alleles reach the assembly or the SNP distances (they must not) |

Reads mimic R10.4.1 SUP reads: about 0.4% substitutions, 0.2% random insertions and deletions, and homopolymer
length errors that grow with the run length; log-normal lengths (N50 about 8 kb); 15% of the reads come from
an unrelated sequence, which baiting must remove (in `numt`, the nuclear reads instead). Real Nanopore errors
are more systematic than these; the real-data comparisons in the records complement the simulations.

In `heteroplasmy` and `numt` each read name ends with its origin (`_major`, `_minor` or `_offtarget`;
`_organelle`, `_nupt<insert>` for a nuclear read covering at least 100 bp of an insert, or `_nuclear`), which
`score.py` counts among the reads BACoN baited (`1_extracted/`) and filtered (`2_filtered/`). The truth genome
of a heteroplasmy sample is its major haplotype (the minor one is in `truth/minor/`); `truth.json` lists each
sample's minor fraction and minor sites, and, for `numt`, the copied segments and every insert with its
substitutions in plastome coordinates (the nuclear sequences are in `truth/nuclear/`). These two scenarios
have their own depths (`--depth` applies to the other three) and their own random streams: adding them left
the files of `plastid`, `linear` and `rdna` byte-identical.

SNP distances are compared with the truth under three conventions: `events` (each mutation once; a SNP in an
inverted repeat or in every copy of a tandem unit is one event), `sites` (each changed position) and
`no_repeat` (SNPs outside repeats only). The convention that fits a method best says how it treats repeats.

## Running

```bash
conda activate BACoN
bash validation/run_validation.sh /tmp/bacon_validation 16          # plastid, linear, rdna
bash validation/run_validation.sh /tmp/bacon_validation_new 16 new  # heteroplasmy, numt
```

Each takes a few minutes (the second about 10: Flye on the mixtures, and generating about 170 Mb of reads).
`python validation/simulate.py FOLDER --scenario new` writes only the new scenarios; `all` (the default) is
still the original three, so earlier records reproduce unchanged. To add a record, copy `versions.txt` and `SCORES_*` into a new
`results/<date>_<version>/` folder with a `SUMMARY.md`. Records are never rewritten: an error found later is
noted as an erratum in a new record.
