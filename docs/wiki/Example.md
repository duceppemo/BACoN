# Example

A small simulated dataset, bundled with BACoN, to check an installation and to see every output in a few
seconds. [`example/make_example.py`](https://github.com/duceppemo/BACoN/blob/main/example/make_example.py)
generates it (standard library only; the same files every time):

- a 30 kb circular reference ("organelle");
- four samples of about 50x of simulated Nanopore reads (1–20 kb, 1.5% errors), each mixed with 20% of reads
  from an unrelated sequence that baiting must remove;
- known SNPs: **alpha** is identical to the reference, **beta** has 6 SNPs, **gamma** has beta's 6 and 4 more,
  **delta** has 10 of its own.

| | Reference | alpha | beta | gamma | delta |
|---|---|---|---|---|---|
| **Reference** | 0 | 0 | 6 | 10 | 10 |
| **alpha** | 0 | 0 | 6 | 10 | 10 |
| **beta** | 6 | 6 | 0 | 4 | 16 |
| **gamma** | 10 | 10 | 4 | 0 | 20 |
| **delta** | 10 | 10 | 16 | 20 | 0 |

## Running it

```bash
conda activate BACoN
bash example/run_example.sh example_output      # optional: number of threads as a second argument
```

The script generates the data in `example_output/data/`, runs BACoN with the default settings into
`example_output/bacon/`, and compares the SNP distances with the table above:

```
OK: all 6 pairwise SNP distances match the truth (example_output/bacon/4_compared/ska/snp_distances.tsv)
```

## The report

Open `example_output/bacon/report.html` in a browser. The report of this example is also online:
[view it](https://duceppemo.github.io/BACoN/reports/example_report.html)
or [download it](https://raw.githubusercontent.com/duceppemo/BACoN/main/docs/reports/example_report.html).
(Local paths were removed from it before publication: the program paths start with `$CONDA_PREFIX`.)

**Overview and samples.** All four samples assembled: about 250 of the 300 reads of each sample matched the
reference (the others are the off-target reads), and after filtering the depth is about 48x (Figure 1). The
templated assembly is 30,000 bp for each sample, as long as the reference, with no `N` base (Figure 2).

![Overview and samples](https://raw.githubusercontent.com/duceppemo/BACoN/main/docs/images/example_report_overview.png)

**Tree.** alpha sits with the reference; beta and gamma form a clade (their 6 shared SNPs); delta is on its
own branch.

![Tree](https://raw.githubusercontent.com/duceppemo/BACoN/main/docs/images/example_report_tree.png)

**SNP distances.** The heatmap, in tree order, shows the pairwise distances; the reference and alpha form a
group of identical genomes (the coloured band).

![SNP distances](https://raw.githubusercontent.com/duceppemo/BACoN/main/docs/images/example_report_distances.png)

**Genome map and VCF.** `example_output/bacon/4_compared/ska/snps.vcf` lists the 20 planted SNPs at their
positions in the reference, with the genotype of each sample; the genome map shows where they are. Every
genome has a call at each of them, and no assembly has an `N` base, so there is no `N` track.

![Genome map](https://raw.githubusercontent.com/duceppemo/BACoN/main/docs/images/example_report_map.png)

**Methods and run.** The methods paragraph describes what was run, with the program versions; the run section
records the command, the reference's MD5 and every program used.

![Methods and run](https://raw.githubusercontent.com/duceppemo/BACoN/main/docs/images/example_report_methods.png)

## Adding sample metadata

A metadata file adds columns to the samples table and colours the tree and the heatmap by one of them
([Usage](Usage#inputs)). Save this as `metadata.csv` (a TSV works too):

```
sample,group,year,comment
alpha,A,2021,Identical to the reference by design
beta,A,2021,Carries six planted SNPs
gamma,B,2022,"Shares beta's six SNPs, plus four of its own"
delta,NA,2023,Ten SNPs of its own; group unknown
```

and run the same command again with `--metadata metadata.csv`: nothing is redone, and the report is rebuilt.
`group` colours the figures (the first column that can be: `year` could too, with `--color-by year`; `comment` is
free text). In the tree, each leaf gets a circle in the colour of its group and the group after its name; delta,
without a value, gets a hollow circle. In the heatmap, a second band gives the group of each genome, and under
the list of identical genomes a table counts, for each group of identical genomes, the genomes of each value.
The groups of identical genomes are then drawn in grey, so that colour always means the metadata.

For real data, see the [Tutorial](Tutorial) (28 potato cultivars) and
[its report](https://duceppemo.github.io/BACoN/reports/tutorial_potato_report.html).
