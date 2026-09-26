# BACoN

**B**ait, **A**ssemble and **Co**mpare **N**anopore reads. BACoN takes Nanopore reads from many samples and a
reference sequence, typically an organelle genome (chloroplast, mitochondrion) or another region of interest,
and for each sample:

1. **baits** the reads that match the reference (minimap2, or BBDuk k-mers);
2. **filters** them by length and quality, and caps the depth (Filtlong);
3. **assembles** them: a templated consensus of the reads aligned to the reference (samtools, the default), or a
   de novo assembly (Flye or myloasm);

then **compares** the samples: SNPs (SKA2 split k-mers by default, or a Parsnp core-genome alignment), a matrix
of pairwise SNP distances and a tree (FastTree or IQ-TREE). Everything is summarized in an HTML report, and in
MultiQC sections.

BACoN was designed for genome skimming: low-coverage sequencing of total DNA, where the organelle reads are
abundant enough to assemble. Enriching the target (targeted sequencing, adaptive sampling) gives more depth
and better assemblies.

| Page | Contents |
|---|---|
| [Installation](Installation) | conda environment, pip, checking the installation |
| [Usage](Usage) | inputs, all options, resuming, performance |
| [Methods](Methods) | what each step does, choosing an assembler and a SNP method, limits |
| [Tutorial](Tutorial) | a complete analysis of 28 potato cultivars from public data |
| [Outputs](Outputs) | every file and column |
| [Validation](Validation) | how the tools were chosen, with simulated data of known truth and real data |
| [FAQ](FAQ) | troubleshooting |
| [Development](Development) | tests, continuous integration, releases |
