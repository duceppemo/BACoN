# FAQ

**A sample failed. Where do I look?**
`summary.tsv` gives the step (`failed (assemble)`) and the reason in `Note`; `logs/<step>/<sample>.log` has
the commands and the programs' messages. A failed sample does not stop the others, and is retried when BACoN
runs again in the same output folder.

**"no reads matched the reference"**
The sample has no read aligning to the reference: wrong reference, a failed library, or a barcode with no
reads. With a distant reference, try `-b bbduk`.

**"No SNP site is shared by all the genomes: no tree"**
With `--snp-method ska`, core SNPs must have their context in every genome. One fragmented or incomplete
assembly is enough to remove every site. Look at `Contigs` and `Length_vs_reference` in `summary.tsv`, remove
the poor samples, or use `--ska-min-freq 0.9` to keep SNPs missing from a few genomes.

**"No tree: no genome differs from the others at any position compared"**
The VCF of the SNPs (all of them, whatever `--ska-min-freq`) is empty: the genomes are identical at every
position SKA2 compared, so there is nothing to build a tree on, and `--ska-min-freq` would not change that. The
distances are all 0.

**What does `--ska-min-freq` below 1 do?**
SKA2 then keeps the SNP sites that some genomes lack (they get a gap in the alignment), as long as at least that
fraction of the genomes has them: a pan-genome SNP alignment instead of core SNPs. It helps when core SNPs are few
or none: fragmented assemblies, gaps, indels or other SNPs within 15 bases of a site, genomes of another
structure. The cost is missing data: each pair of genomes is compared on the sites both have, and a genome
missing many sites looks closer to the others than it is. `snps.vcf` is the same whatever the value. Details in
[Methods](Methods#4-comparison).

**The de novo assembly is 1.3–2 times the reference length.**
Reads from other genomes were baited with the target: mitochondrial and nuclear sequences similar to the
plastome (for example plastid DNA inserted in the mitochondrion), or contamination. The extra contigs are
usually at low depth. The templated assembly is not affected; for de novo assemblies, check the assembly
graph in `3_assembled/assembly_graphs/`.

**My chloroplast assembly is in three contigs, not one circle.**
With reads shorter than the inverted repeat (about 25 kb), the assembler cannot tell how the single-copy regions
connect through the repeat: the large single copy, the repeat (once) and the small single copy come out
separately. It does not affect SNP calling with SKA2. Longer reads give a circular assembly.

**Flye made a small circular contig of my linear target.**
Flye collapses tandem repeats (such as rDNA arrays) into one circular unit, dropping the flanks
([Validation](Validation)). Use the templated assembly, or `-a myloasm`, which keeps them linear.

**Why are there `N` in the templated assembly?**
Where fewer than three reads cover the reference, or the reads disagree on a base (often inside an insertion), samtools
writes `N`. `N_bases` in `summary.tsv` counts them. A de novo assembly resolves insertions.

**Which reads can I use?**
Any Nanopore reads, fastq or fasta. Dorado or Guppy 5+ SUP/HAC reads work best. Flye's default `nano-hq` mode
is meant for R10/Q20 reads (<5% error); for R9 reads basecalled with Guppy < 5 (10–20% error), add
`--read-type nano-raw`. The templated assembly works with either. The tutorial's older MinION reads (mean
Q14–16) worked with the defaults, including `-a flye`.

**Can I compare samples sequenced with different chemistries?**
Yes, but systematic consensus errors differ between chemistries and add to the distances between samples
from different runs. The templated assembly is the least affected.

**How do I compare with a public genome?**
Give its fasta file with `--add-genomes`: finished genomes (published plastomes, earlier assemblies) join
the comparison as they are, named after their file, without baiting or assembly.

```bash
bacon -r ref.fasta -i reads/ -o out/ --add-genomes published/*.fasta
```

**How do I resume, or change one parameter?**
Run the same command again; see [Usage](Usage#resuming-and-changing-parameters).

**Where did Porechop, Shasta, Rebaler, Snippy and PhaME go?**
They were removed in 0.3, with the evidence in [Methods](Methods#choices-that-changed-in-03) and
[Validation](Validation). BACoN 0.2 is still available from the
[releases](https://github.com/duceppemo/BACoN/releases).
