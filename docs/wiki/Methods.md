# Methods

## 1. Baiting

The reads of each sample are compared with the reference and those that match are kept
(`1_extracted/`).

- **minimap2** (default): `minimap2 -x map-ont --secondary=no`; a read is kept when it has any alignment to
  the reference. The read is copied unchanged (not trimmed to the aligned part), with its qualities.
- **BBDuk** (`-b bbduk`): a read is kept when it shares one k-mer (k = 31 by default, one mismatch allowed: `--hdist`) with
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
| Output | one sequence per reference sequence with reads (one without any read is left out, with a note); positions covered by fewer than three reads are `N` | contigs; circular ones flagged | contigs; circular ones flagged |
| Accuracy of the consensus | highest | high; some systematic errors on real data | high |
| Structure (circularity, rearrangements, inverted repeats) | not shown: follows the reference | shown | shown |
| Insertions absent from the reference | small ones (up to a few hundred bp) | all | all |
| Tandem repeats (e.g. rDNA arrays) | follows the reference | may collapse into one circular unit | kept linear |
| Needs | a reference close to the samples (same species) | nothing but depth (about 20x or more) | depth |
| Assembly time, 28 potato plastomes (`-t 48 -p 12`) | 13 s | 10 min | 35 s |

**Which one?** For comparing samples of one species against a good reference, which is the usual aim of
BACoN, use the templated assembly (the default): in every test it gave SNP distances as good as or better than
the de novo assemblies, and it is the fastest ([Validation](Validation)). Use a de novo assembler to look at
the structure of the genomes, when the reference is distant, or to check the templated results
independently: on the tutorial data, Flye and the templated assembly gave identical SNP distances for all
406 pairs of genomes (the 28 cultivars and the reference). Between the two de novo assemblers, Flye gave more complete plastomes on real data;
myloasm is faster and does not circularize linear tandem arrays.

Limits of the templated assembly:

- It shows the sample only where its reads align to the reference; a region absent from the reference
  cannot appear, and large insertions are left out. Insertions and deletions of up to a few hundred bases
  are kept, such as the 241 bp that distinguishes potato cytoplasm types ([Tutorial](Tutorial)), but the
  bases of an insertion the reads do not agree on are `N`: Flye resolved that insertion fully.
- Where fewer than three reads align, or the reads disagree, the consensus has `N` (`N_bases`, and a note in
  `summary.tsv`).
  This is `--template-gaps n`, the default; with `--template-gaps reference`, uncovered positions at the ends of
  each sequence get the reference's bases instead (inner gaps stay `N`).
- In an inverted repeat both copies receive the same reads; a difference between the two copies of one
  sample cannot be seen (the two copies are usually identical in chloroplasts).

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
  their first 30 bases so that SNPs next to the start of the sequence are not lost, and so is the reference
  and the genomes of `--add-genomes`
  when more than half of the assemblies of the samples have a circular contig (Flye and myloasm report them; templated
  assemblies follow the reference and are never flagged circular, so with them the SNPs within 15 bases of the
  ends of the reference are not compared: a split k-mer needs 15 bases on either side of its middle base).
- **Parsnp** (`--snp-method parsnp`) aligns the core genome of the assemblies to the reference. BACoN runs it
  with `-c`, which turns off Parsnp's divergence (MUMi) filter: by default Parsnp drops genomes it finds too
  divergent, which in the validation removed a low-depth sample. A genome much shorter than the reference can
  still make Parsnp fail, or be kept and reduce the core genome of every genome to what it covers: BACoN warns
  when the core genome is less than half of the reference. SNPs inside inverted repeats are missed, and SNPs
  next to large deletions may be. SNPs within about 30 bases of the ends of the reference are not compared
  either: Parsnp's locally collinear blocks end at their last maximal unique match. HarvestTools writes to
  `snps.vcf` the SNPs Parsnp filtered out of its alignment too, with a `FILTER` other than `PASS` (`ALN`: more
  than 20 indels in the 100 bp window around the SNP; `CID`: that window less than 50% identical; `LCB`: a block
  shorter than 200 bp; `IND`, `N`: a column with an indel or an N; several joined by `:`, `CID:ALN`). Those
  SNPs are not in the SNP alignment, so not in the distances or the tree; the report shows the `PASS` records
  only and says how many others `snps.vcf` has. SKA2's VCF has no filter (`.`): every record is shown.

**Core or pan-genome SNPs (`--ska-min-freq`, SKA2 only).** SKA2 keeps a SNP site when at least this fraction of
the genomes (the reference and the genomes of `--add-genomes` included) has its split k-mer: the same 15 bases on
either side of it. At 1 (default), every genome must have it: a core-SNP alignment, with a nucleotide for every
genome in every column. Below 1 (0.5: at least half of the genomes), the sites missing from some genomes are kept,
with a gap (`-`) for those genomes: a pan-genome SNP alignment, the role kSNP3 had in BACoN 0.1. The results go
to `4_compared/ska_<value>/` instead of `4_compared/ska/`.

A genome lacks the split k-mer of a site when:

- the site is within 15 bases of the end of a contig or of a run of `N`: fragmented de novo assemblies (Flye or
  myloasm at low depth), gaps of a templated assembly;
- another SNP or an indel lies within 15 bases of it in that genome (its flanks differ, so its split k-mer is
  another one);
- the region is missing or rearranged in that genome (a deletion, a published genome with another structure).

With core SNPs, one such genome removes the site for all; when that happens to every site, no tree is built
("No SNP site is shared by all the genomes"). A lower `--ska-min-freq` keeps those sites, at a cost: the
alignment has missing data, and as the distances count only the positions where both genomes have a nucleotide,
each pair is compared on its own set of sites. A genome missing many sites then looks closer to the others than
it is, and distances between different pairs are less directly comparable. FastTree and IQ-TREE treat the gaps
as missing data. In the [validation](Validation), `--ska-min-freq 0.5` gave the same exact distances as core SNPs
on the simulated plastid set.

`--ska-min-freq` does not change `snps.vcf`, which lists every SNP relative to the reference (from `ska map`),
nor genomes that do not differ at any position compared (an empty VCF): those give no tree whatever the value.
Parsnp has no equivalent: it aligns the core genome.

**Distances.** `snp_distances.tsv` counts, for each pair of genomes, the positions of the SNP alignment where
both have a nucleotide and they differ.

**Tree.** FastTree (GTR, SH-like local supports from 100 resamples; default) or IQ-TREE (`--tree iqtree`: ModelFinder and 1000
ultrafast bootstraps; with fewer than four distinct sequences IQ-TREE makes no bootstrap, and its maximum-likelihood tree
is shown without supports), built on the SKA2 SNP alignment or the Parsnp core-genome alignment. The tree is rooted
at the midpoint of its longest path, ladderized, and drawn as SVG. With SKA2, branch lengths are in
substitutions per variable site.

At least three assemblies are needed for a comparison. When no SNP site is shared by all the genomes (for
example with fragmented assemblies), or when the genomes do not differ at any position compared, the distances are
written but no tree is built.

## 5. SNP effects

With an annotation of the reference (`--annotation`, or a GenBank reference), the report places each SNP of
`snps.vcf` on the reference's genes. BACoN reads GenBank feature tables (`gene`, `CDS`, `tRNA`, `rRNA`, their
locations with `join`, `complement`, `order`, partial ends and trans-splicing, `/codon_start`,
`/transl_table`, `/transl_except`, `/pseudo`, quoted values spanning lines with `""` for a quote) and GFF3
(features linked by `ID` and `Parent`, the phase column, `pseudogene` features and the `pseudo`, `gene_biotype`
and `biotype` attributes, `transl_except`); a CDS or RNA feature belongs to the gene feature of the same name
that overlaps it, or is a gene of its own. The exons of a gene are the parts of its CDS and RNA features; a GFF3
RNA feature with `exon` features of its own (NCBI writes a spliced tRNA, such as plastid *trnK*, *trnL*,
*trnV*, *trnI*, *trnA* or *trnG*, as one line from its first base to its last) has those exons only, so its
intron is an intron in either format (*matK*, inside the intron of *trnK*, is `CDS / intron`). A gene is named after its `gene` qualifier, else `locus_tag`,
`product`, `Name`, `gene_id` or `ID`. Older RefSeq records (plastomes among them) have genes without a `/gene`
symbol, only a `/locus_tag` and a `/product` (recent records usually have the symbol): such a gene keeps its
locus tag as its identifier, but the report shows it as `locus_tag (product)`, `LK299_pgr007 (23S ribosomal
RNA)` in NC_008096.2 (potato), in the SNP table's Gene column, the hovers, the intergenic contexts (`intergenic
between X and Y`) and the summary of the genes with the most SNPs (`LK299_pgr007 (23S ribosomal RNA; 19)`); on
the map, where labels must be short, by its product when that is at most 12 characters (`tRNA-Val` for
`OrsajCt141` of NC_001320.1, rice), else by its identifier (`LK299_pgr007`). A gene has a symbol when a feature
of its has `/gene`; in GFF3, `gene=`, or on the gene feature a `Name` that is not an identifier: its locus tag,
its product, its `gene_id` or its ID, or its ID without the type prefix (`gene-`, `gene:`) when the gene has a
locus tag or the `Name` looks like one (`PREFIX_number`; NCBI names a gene without a symbol after its locus tag,
whereas a third party's `ID=gene-matK;Name=matK` keeps *matK* as its symbol; the `Name` of a CDS, a protein
accession, is never a symbol). No symbol is made up from a product. Ensembl's `ncRNA_gene` features are genes
(their tRNA and rRNA transcripts are named after them), and a CDS or RNA feature without a `Parent` belongs to
the gene feature on its strand that contains it and shares its locus tag or symbol (or, when it has neither, to
any gene feature containing it: NCBI writes the rRNA of a gene with `gene_biotype=other` without a `Parent`).
The lines of a GenBank value continued on the next line are joined with a space, except after a hyphen or
after a comma followed by a digit, where the flat file broke a word (`2,` / `6-diaminopimelate`). The parts of
a CDS are read in the order of translation: as the GenBank location lists them,
and for the lines of a GFF3 CDS (sharing an `ID`, or without one, a `Parent`) by their `part=` numbers (NCBI)
when given, else in file order (NCBI lists them 5′ to 3′), except for lines listed by ascending coordinate
(Ensembl, a sorted file): those of a −-strand CDS are read in descending order, and those of a CDS across the
origin (spanning more than half the sequence, from its first base to its last) start at the high part on the +
strand and at the low part, read downwards, on the − strand. The sequence length comes from `##sequence-region`,
else from the reference; without either (an annotated sequence matched to the reference by nothing but being
the only one), a CDS across the origin listed by coordinate is read in coordinate order. A single line ending
beyond the length of a sequence flagged `Is_circular=true` (Bakta writes a feature across the origin so, with
end = its end + the length) is split into its two parts, and a line lying wholly beyond that length (NCBI
writes the 5′ exon of the *Epifagus* NC_001568.1 *rps12* at 82,777–82,890 on its 70,028 bp sequence:
12,749–12,862) is moved back by the length. A trans-spliced CDS (`exception=trans-splicing`, or parts on both
strands) whose lines are listed by ascending coordinate without `part=` numbers has no known order unless the
file lists the parts of its other −-strand features 5′ to 3′ (NCBI); in a sorted file (Ensembl) it gives no
effects, with a warning, rather than be read in a guessed order (descending, for potato's IRb *rps12*
LK299_pgp043, would read exon 3 first). So the trans-spliced plastid *rps12*, whose 5′ exon
lies in the LSC on the other strand, and
a gene across the origin of a circular sequence are translated correctly from either format; a CDS or RNA
feature across the origin (or trans-spliced) without a gene feature is a gene made of its parts, not of the
whole sequence between them. A CDS with a part on another sequence (`join(X12345.1:1..100,201..500)`) keeps
its parts on this one, but its frame there is unknown: it gives no effects, with a warning.

For a SNP inside a coding sequence, BACoN rebuilds the coding sequence from the reference (the exons in the
order of translation, reverse-complemented on the − strand, from the base given by `codon_start`), finds the
codon containing the SNP, substitutes the alternate allele (complemented on the − strand) and translates both
codons with the CDS's translation table (`/transl_table`, otherwise table 11, the bacterial and plastid code,
or table 1 for a sequence said to be a mitochondrion by its GenBank `/organelle` or NCBI GFF3 `genome=`: plant
mitochondrial records give no `/transl_table`, and the INSDC default is table 1; tables 1, 2, 3, 4, 5, 9, 11, 13
and 14 are known, as [NCBI defines them](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi), table 1
being the standard code with ATG as its only start codon; others, such as 16, 21, 22 or 23, use the standard
code with table 11's start codons, with a warning). The effect is `synonymous`, `missense`,
`nonsense` (a stop codon gained), `stop lost`, `stop retained` (a stop codon changed into another), and for the
initiation codon of a complete CDS `start lost` or `start retained` (the new codon is another start codon of
the table, such as GTG in table 11). No effect is given when the codon is incomplete (a partial CDS), contains
`N`, or when the VCF's reference allele does not match the reference, nor for a codon with a translational
exception (`/transl_except`: a selenocysteine or pyrrolysine codon, an edited codon, a stop codon completed by
polyadenylation, a codon split by an intron, `pos:join(...)`), which is not a stop codon although the table
says so: the SNP is still reported in the CDS, without an effect. RNA editing (plastid ACG start codons, for example) is not modelled. Pseudogenes get no
effect. Each alternate allele, and each of two overlapping coding sequences, gets an effect of its own, but the
same effect through two coding sequences of a gene (the 5′ exon shared by the two *rps12* of a plastome, the
two products of a ribosomal frameshift) is given once. A base read twice through a −1 ribosomal frameshift,
written as overlapping parts (`join(66..327,327..1228)` in the F plasmid NC_002483.1, *dnaX* of *E. coli*),
is changed in every copy: one effect per codon it changes, each codon with both copies changed.

Outside coding sequences, the context is the gene's type (tRNA, rRNA), `intron` when the position lies
between two exons of a gene, or `intergenic between X and Y`, the nearest genes on either side (around the
origin when the annotation says the sequence is circular). A trans-spliced gene (plastid *rps12*) has no
intron between its distant parts.

The LSC/IRb/SSC/IRa band of a plastome is derived from the annotated inverted repeats (`repeat_region` with
`/rpt_type=inverted`, or any region feature whose note starts by naming an inverted repeat, IRa or IRb, such
as `IRb`, `inverted repeat A` or `inverted repeat region IRa`; a note merely mentioning a copy, `... in IRA`
or `junction LSC-IRB`, or mentioning a junction or a border, `IRB/SSC junction` or `IRA-SSC border`, does not
annotate it, nor does a feature shorter than 50 bp, so the junctions potato marks on 2 bp do not grow the
copies), when there are two of at least 500 bp: the single-copy regions
are the gaps between them, the larger one being the LSC, and the repeat following the LSC is IRb unless the
annotation names them; one region may span the origin, a repeat too (`join(x..length,1..y)`, or two features
named alike meeting at the origin). A repeat annotated twice (a `repeat_region` and a `misc_feature` over the
same stretch) counts once; with more than two repeats annotated, BACoN takes the two `repeat_region` with
`/rpt_type=inverted`, else the pair named IRa and IRb, else the two longest.

Most older RefSeq plastomes do not annotate their inverted repeats, so BACoN also looks for the large inverted
repeat in each reference sequence itself, with or without an annotation (standard library only, linear in the
length: about 0.1 s for a plastome, 1 s for 2 Mb; sequences longer than 2 Mb, bacterial chromosomes, are not
searched). Every fourth 32-mer of the sequence is indexed (a 32-mer at several indexed positions, or with a
base other than ACGT, is left out), every 32-mer of the reverse complement is looked up, and each match is a
seed (*i*, *j*): the 32-mer at *j* read on the other strand is the one at *i*. The seeds of an inverted repeat
share an antidiagonal (*i* + *j* constant; an indel between the copies shifts it by its length), so the
antidiagonal with the most seeds and those within 500 of it are chained in order (*i* increasing, *j*
decreasing; a gap of more than 2 kb without a seed ends the chain; a seed on a parallel antidiagonal is only
taken when the current one has no seed within 2 kb after it, which keeps a short duplication inside the repeat
from being read as an indel), and the chain goes on past a larger indel, of up to 3 kb in one copy, when the
seeds resume beyond it on another antidiagonal (an insertion of 600 bp or of 2.5 kb in one copy does not cut
the repeat; a longer one, or a run of N longer than 2 kb, does, and the longer part is the repeat). The
differences between the copies are counted between the seeds: mismatches, indels and runs of bases other than
ACGT, each counted once whatever its length (an insertion, or a scaffold gap of 300 N, is one difference). The
stretches between two seeds are compared base by base when they have the same length, and aligned (an edit
distance within a band of their length difference plus 16 bases) when their lengths differ or when the
comparison base by base finds more than three mismatches: two compensating indels, a base lost and another
gained 300 bp further, look like 225 mismatches base by base but are two differences (NC_007144.1, cucumber,
whose copies are 99.9% identical, has such a pair). The chain is then trimmed to its part scoring most (a base
+1, a difference −20), which drops a diverged flank (the Arabidopsis mitochondrion NC_037304.1 has a 6,590 bp
repeat, 100% identical, inside flanks far less so), and extended base by base at both ends, through a mismatch
when the 12 bases after it match, so the copies' coordinates are those of the maximal match (they can differ
by one to eight bases from a published junction where the flanking bases happen to match; a copy may end across
the origin, as the tomato NC_007898.3 IRa does by one base). A repeat is accepted when each copy is at least
5,000 bp (plastid IRs are 10–30 kb; a lineage that lost one copy, or the 3 kb repeat of the bundled example,
gets no band from the sequence), the copies are at least 99% identical and do not overlap; copies that abut (the
two halves of a palindrome; the *Toxoplasma* apicoplast NC_001799.1 has its copies abutting across the origin)
are a repeat without a region between them, so without a band. A repeat across the origin of the circular
sequence, or ending exactly at it with matching 32-mers beyond, is searched again from a point between the
copies, and that result is taken when it is the same pair of copies. The regions are then derived as from
annotated repeats, named by convention (IRb follows the LSC).

The band is drawn only when the four regions have the quadripartite layout of a plastome, whether the repeats
are annotated or detected: the two repeats make at least 5% of the sequence, the larger single-copy region is at
most 200 kb (among 64 plastid genomes tested, from the 11 kb *Pilostyles* to the 218 kb *Pelargonium* plastome,
the 48 with an inverted repeat of 5 kb have 10–70% of their sequence in the repeats and a LSC of at most 136 kb;
the inverted rRNA operons of *Helicobacter pylori* NC_000915.1 are 1.3% of 1.67 Mb, and the 16.9 kb repeat of the
maize mitochondrion NC_007982.1 5.9% of 570 kb with 481 kb between the copies), and the record is not said to be
something else: a GenBank `source` feature whose `/organelle` is not a plastid (`mitochondrion`), or an NCBI GFF3
`region` whose `genome=` is not one (`chromosome`, `mitochondrion`), gets no plastome names. Such a repeat,
annotated or detected, is recorded in `run_info.json` without regions and with a note saying why, and the log
says so.

The annotated repeats keep priority: when the annotation marks two inverted repeats of at least 500 bp that make
a plastome layout, the band follows them (their coordinates) as long as each copy of the repeat found in the
sequence lies at least 80% inside an annotated copy (a detection cut short by a large insertion or a run of N
confirms the annotation). The sequence wins, with a warning, when the detected copies lie elsewhere, or when the
annotated repeats are not a plastome layout. NC_001879.2 (tobacco) annotates its LSC as `inverted repeat B`,
its SSC as `inverted repeat A` and its IRa as a third feature: its features give no layout at all, which the
report says, and the band follows the sequence. The map's caption, the Methods paragraph (when the report shows
a band) and the log say which source was used, with the detected copies' size and identity (`the inverted
repeat detected in the reference sequence (two copies of 25,593 bp, 100% identical)`), and `run_info.json`
records, for each reference sequence, the regions, their source and the detected copies, or `"none"`. On
NC_008096.2 (potato) the detected copies are exactly the annotated IRb 85,738–111,330 and IRa 129,704–155,296;
on NC_000932.1 (*Arabidopsis*) 26,264 bp, the published size; on NC_007144.1 (cucumber) 25,193 and 25,189 bp,
the copies minimap2 aligns.

## Choices that changed in 0.3

| 0.2 | 0.3 | Why |
|---|---|---|
| Porechop | removed | unmaintained, slow; reads are trimmed by the basecaller |
| Shasta | removed | on real data, some assemblies collapsed to a fraction of the genome; on simulations it left overlapping ends and lost part of an inverted repeat |
| Rebaler | samtools consensus | Rebaler is unmaintained (2019); samtools consensus was far more accurate in every test |
| Snippy | SKA2 | Snippy (2020) can no longer be installed with current assemblers; it misses SNPs in inverted repeats and was the noisiest on real data |
| PhaME | removed | its dependency check compares versions as decimals and rejects samtools 1.10 to 1.29 (1.21 reads as older than 1.3; observed with 1.16 and 1.21) |
| Parsnp without `-c` | Parsnp with `-c` | no longer drops divergent genomes (a genome much shorter than the reference can still make Parsnp fail, or shrink the core genome; BACoN warns) |
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
