# Outputs

```
OUTPUT/
├── summary.tsv                 one line per sample: reads, depth, assembly, status, notes
├── run_info.json               version, command, settings, program versions, samples, comparison
├── bacon.log                   the log of every run in this folder
├── reference.fasta             the reference used (uncompressed copy)
├── 1_extracted/<sample>.fastq.gz          baited reads (and <sample>.bam with --keep-bam)
├── 2_filtered/<sample>.fastq.gz           filtered reads
├── 3_assembled/
│   ├── all_assemblies/<sample>.fasta      the assembly of each sample
│   ├── assembly_graphs/<sample>.gfa/.png  assembly graphs (de novo assemblers; .png with Bandage)
│   └── <sample>/                          the assembler's own folder
├── 4_compared/<method>/
│   ├── snp_distances.tsv       pairwise SNP distances, Reference first
│   ├── tree.nwk                midpoint-rooted tree (Newick)
│   ├── tree.svg                picture of the tree
│   └── ...                     the method's own files (alignments)
└── logs/<step>/<sample>.log    the commands run and the programs' messages
```

## summary.tsv

| Column | Description |
|---|---|
| `Sample` | Sample name |
| `Status` | `ok`, or `failed (step)` with the reason in `Note` |
| `Raw_reads`, `Raw_bases` | Input reads and bases (`NA` with `-b bbduk` when BBDuk does not report them) |
| `Baited_reads`, `Baited_bases`, `Baited_pct` | Reads matching the reference, and their share of the input bases (%) |
| `Filtered_reads`, `Filtered_bases`, `Filtered_N50` | Reads kept by Filtlong |
| `Est_depth` | Filtered bases / reference length |
| `Contigs`, `Assembly_length`, `Largest_contig`, `Assembly_N50` | Assembly statistics |
| `Circular_contigs` | Contigs the assembler reports as circular (Flye, myloasm; `NA` for the templated assembly) |
| `Length_vs_reference` | Assembly length / reference length |
| `Assembly_depth` | Mean depth reported by Flye |
| `N_bases` | Templated assembly: bases called `N`, where no read covers the reference or the reads disagree (often inside insertions) |
| `Note` | Warnings: low depth (below 20x), assembly length outside 0.8–1.2 times the reference, `N` bases, not in the tree, or why the sample failed |

## Assemblies

`3_assembled/all_assemblies/<sample>.fasta` holds one sequence per contig, named `<sample>_<contig>`. Contigs
that the assembler reports as circular have ` circular=true` in their header.

- **Templated** (`samtools`): one sequence per reference sequence, in reference coordinates plus the
  sample's insertions and minus its deletions. Positions that no read covers, and bases the reads do not agree
  on (often inside an insertion), are `N`.
- **De novo** (`flye`, `myloasm`): contigs as assembled. A circular genome can start anywhere and on either
  strand. Chloroplasts often come out as three contigs (large single copy, inverted repeat, small single copy)
  when the reads are shorter than the inverted repeat.

## Comparison

`snp_distances.tsv` is a square matrix of the number of positions where two genomes have different
nucleotides (A, C, G, T; gaps and N are ignored), with the reference as `Reference`.

- **SKA2** (`4_compared/ska/`, or `ska_<min-freq>/`): `ska.snps.fasta` is the alignment of the variable sites.
  A SNP in an inverted repeat is counted once (its two copies share the same split k-mer).
- **Parsnp** (`4_compared/parsnp/`): `parsnp.core.fasta` is the core-genome alignment, used for the tree;
  `parsnp.snps.fasta` the SNP sites; plus Parsnp's own files.

`tree.nwk` is rooted at the midpoint of the longest path and ladderized; internal labels are the supports
(SH-like local supports for FastTree, ultrafast bootstraps for IQ-TREE). `tree.svg` draws it with a scale in
substitutions per site. When no SNP site is shared by all the genomes, the distances are written but no tree
is built; lower `--ska-min-freq` (see [FAQ](FAQ)).

## run_info.json

The provenance of the last run: BACoN version, command line, start time and duration, all settings, the path
and version of every program used, each sample's files and status, and the comparison's result files.
