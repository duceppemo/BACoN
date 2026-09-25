# Methods

## 1. Baiting

The reads of each sample are compared with the reference and those that match are kept
(`1_extracted/`).

- **minimap2** (default): `minimap2 -x map-ont --secondary=no`; a read is kept when it has any alignment to
  the reference. The read is copied unchanged (not trimmed to the aligned part), with its qualities.
- **BBDuk** (`-b bbduk`): a read is kept when it shares one k-mer (k = 31 by default, up to two mismatches) with
  the reference. Useful when the reference is only distantly related to the samples.

Baiting reads every input read once; everything after works on the baited reads only. The share of bases
baited (`Baited_pct`) is a measure of enrichment: a few percent in total DNA of plants, much more after
targeted sequencing.

## 2. Filtering

Filtlong discards reads shorter than `--min-read-length` (500 bp), keeps the best `--keep-percent` (95%) of
the bases by length and quality, and caps the depth at `--target-depth` (100x the reference length). Capping
keeps assemblies fast on highly enriched samples without losing accuracy.

Filtlong needs base qualities: reads given as fasta are filtered by BACoN instead, keeping the longest reads of
at least `--min-read-length` up to the target depth.

Adapters are not trimmed: Dorado and Guppy already trim them, and Porechop, used by BACoN 0.2, is no longer
maintained.

## 3. Assembly

| | Templated: `samtools` (default) | De novo: `flye` | De novo: `myloasm` |
|---|---|---|---|
| How | reads aligned to the reference (minimap2), consensus with `samtools consensus -X r10.4_sup` | Flye 2.9 | myloasm |
| Output | one sequence per reference sequence; positions without reads are `N` | contigs; circular ones flagged | contigs; circular ones flagged |
| Accuracy of the consensus | highest | high; some systematic errors on real data | high |
| Structure (circularity, rearrangements, inverted repeats) | not shown: follows the reference | shown | shown |
| Insertions absent from the reference | small ones (up to a few hundred bp) | all | all |
| Tandem repeats (e.g. rDNA arrays) | follows the reference | may collapse into one circular unit | kept linear |
| Needs | a reference close to the samples (same species) | nothing but depth (about 20x or more) | depth |
| Assembly time, 28 potato plastomes (`-t 48 -p 12`) | 17 s | 10 min | 35 s |

**Which one?** For comparing samples of one species against a good reference, which is the usual aim of
BACoN, use the templated assembly (the default): in every test it gave SNP distances as good as or better than
the de novo assemblies, and it is the fastest ([Validation](Validation)). Use a de novo assembler to look at
the structure of the genomes, when the reference is distant, or to check the templated results
independently: on the tutorial data, Flye and the templated assembly gave identical SNP distances for all
406 pairs of samples. Between the two de novo assemblers, Flye gave more complete plastomes on real data;
myloasm is faster and does not circularize linear tandem arrays.

Limits of the templated assembly:

- It shows the sample only where its reads align to the reference; a region absent from the reference
  cannot appear, and large insertions are left out. Insertions and deletions of up to a few hundred bases
  are kept, such as the 241 bp that distinguishes potato cytoplasm types ([Tutorial](Tutorial)), but the
  bases of an insertion the reads do not agree on are `N`: Flye resolved that insertion fully.
- Where no read aligns, or the reads disagree, the consensus has `N` (`N_bases`, and a note in
  `summary.tsv`).
  `--template-gaps reference` copies the reference into uncovered sequence ends only; inner gaps stay `N`.
- In an inverted repeat both copies receive the same reads; a difference between the two copies of one
  sample cannot be seen (such differences are rare: copies are homogenized in chloroplasts).

Limits of the de novo assemblies:

- Below about 15–20x, assemblies break into pieces and may duplicate regions (Validation, sample at 10x).
- When reads are shorter than the inverted repeat of a chloroplast (about 25 kb), the assembly is usually three
  contigs, with the repeat once.
- With ONT reads, all samples share some systematic consensus errors (homopolymers). They cancel between samples
  but add to the distance to the reference.

## 4. Comparison

**SNPs.**

- **SKA2** (default, `--snp-method ska`) splits every genome into split k-mers (two 15-mers around a middle
  base) and aligns the middle bases of the split k-mers found in the genomes; it needs no reference and no
  alignment. `--ska-min-freq 1` (default) keeps SNPs whose context is in every genome: core SNPs.
  Lower values keep SNPs missing from some genomes (a pan-genome SNP alignment with gaps, like kSNP, which
  BACoN 0.1 used). A SNP in an inverted repeat is counted once. Before SKA2, circular contigs are extended by
  their first 30 bases so that SNPs next to the start of the sequence are not lost.
- **Parsnp** (`--snp-method parsnp`) aligns the core genome of the assemblies to the reference. BACoN runs it
  with `-c` so that it keeps every genome: by default Parsnp drops genomes it finds too divergent, which in
  the validation removed a low-depth sample. SNPs inside inverted repeats are missed, and SNPs next to
  large deletions may be.

**Distances.** `snp_distances.tsv` counts, for each pair of genomes, the positions of the SNP alignment where
both have a nucleotide and they differ.

**Tree.** FastTree (GTR, SH-like local supports; default) or IQ-TREE (`--tree iqtree`: ModelFinder and 1000
ultrafast bootstraps), built on the SKA2 SNP alignment or the Parsnp core-genome alignment. The tree is rooted
at the midpoint of its longest path, ladderized, and drawn as SVG. With SKA2, branch lengths are in
substitutions per variable site.

At least three assemblies are needed for a comparison. When no SNP site is shared by all the genomes (for
example with fragmented assemblies), the distances are written but no tree is built.

## Choices that changed in 0.3

| 0.2 | 0.3 | Why |
|---|---|---|
| Porechop | removed | unmaintained, slow; reads are trimmed by the basecaller |
| Shasta | removed | on real data, some assemblies collapsed to a fraction of the genome; on simulations it left overlapping ends and lost part of an inverted repeat |
| Rebaler | samtools consensus | Rebaler is unmaintained (2019); samtools consensus was far more accurate in every test |
| Snippy | SKA2 | Snippy (2020) can no longer be installed with current assemblers; it misses SNPs in inverted repeats and was the noisiest on real data |
| PhaME | removed | broken with any recent samtools (its version check reads 1.21 as older than 1.3) |
| Parsnp without `-c` | Parsnp with `-c` | keeps every genome |
| RAxML 8 | IQ-TREE | maintained; model selection |
| ete3 PDF | SVG, drawn by BACoN | ete3 needs Qt and a display |

The evidence is in [Validation](Validation).

## References

- minimap2: Li H. (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34:3094–3100.
- samtools consensus: Danecek P. et al. (2021) Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008.
- Filtlong: Wick R. https://github.com/rrwick/Filtlong
- Flye: Kolmogorov M. et al. (2019) Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology* 37:540–546.
- myloasm: Shaw J., Marin M.G., Li H. (2026) High-resolution metagenome assembly for modern long reads with myloasm. *Nature Biotechnology*. https://doi.org/10.1038/s41587-026-03053-z
- SKA2: Derelle R., von Wachsmann J., Mäklin T., Hellewell J., Russell T., Lalvani A., Chindelevitch L., Croucher N.J., Harris S.R., Lees J.A. (2024) Seamless, rapid, and accurate analyses of outbreak genomic data using split k-mer analysis. *Genome Research* 34:1661–1673. https://doi.org/10.1101/gr.279449.124
- Parsnp: Treangen T.J. et al. (2014) The Harvest suite for rapid core-genome alignment and visualization of thousands of intraspecific microbial genomes. *Genome Biology* 15:524; Kille B. et al. (2024) Parsnp 2.0: scalable core-genome alignment for massive microbial datasets. *Bioinformatics* 40(5):btae311. https://doi.org/10.1093/bioinformatics/btae311
- FastTree: Price M.N. et al. (2010) FastTree 2. *PLoS ONE* 5:e9490.
- IQ-TREE: Wong T.K.F. et al. (2026) IQ-TREE 3: phylogenomic inference software using complex evolutionary models. *Molecular Biology and Evolution* 43(5):msag117. https://doi.org/10.1093/molbev/msag117
- BBDuk: Bushnell B. BBTools. https://sourceforge.net/projects/bbmap/
