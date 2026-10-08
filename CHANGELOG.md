# Changelog

## Unreleased

### Added
- The LSC/IRb/SSC/IRa band of the genome map no longer needs annotated inverted repeats: BACoN detects the
  large inverted repeat in each reference sequence (two copies of at least 5 kb, at least 99% identical, in a
  sequence of 2 Mb or less; about 0.1 s for a plastome, 1 s for 2 Mb) and uses it when the annotation has no
  inverted repeats (older RefSeq plastomes often annotate none) or when there is no annotation at all (the gene
  track still needs one). The copies' coordinates are those of the maximal match: on NC_008096.2 (potato)
  exactly the annotated ones, on other records within 1–8 bp of the annotated junctions, on NC_007144.1
  (cucumber) the copies minimap2 aligns. Mismatches, indels and runs of N each count as one difference, so an
  insertion of up to 3 kb in one copy, a scaffold gap or compensating indels (cucumber) do not reject the
  repeat, and a diverged flank is dropped (the Arabidopsis mitochondrion NC_037304.1 has a 6,590 bp repeat,
  100% identical, in diverged flanks). The band is only drawn when the regions have the layout of a plastome:
  the repeats at least 5% of the sequence, the larger single-copy region at most 200 kb, and the record not said
  to be a mitochondrion (GenBank `/organelle`, NCBI GFF3 `genome=`) or a chromosome; the inverted rRNA operons
  of a bacterium (*Helicobacter pylori* NC_000915.1) or the repeat of a plant mitochondrion (maize NC_007982.1)
  get no band. Annotated repeats keep priority when the detected copies lie inside them; the sequence wins,
  with a warning, when the annotation names other stretches, or names repeats that give no layout (NC_001879.2,
  tobacco, annotates its LSC as `inverted repeat B` and its IRa twice), or repeats that are not a plastome
  layout. The caption, the Methods paragraph (when the report shows the band) and the log say where the regions
  come from, with the detected copies' size and identity; `run_info.json` records the regions of each reference
  sequence (`reference.regions`), their source and the detected copies, `"none"`, or a repeat that gives no
  band with a `note` saying why.

### Changed
- Each sample has one colour and shape throughout the report, from a single source (`sample_mark`) shared by
  every figure and table. The bar charts of depth and of `N` bases, which were all one colour, now colour each
  sample's bar like the rest of the report and draw its marker before its name: with a colour column, the
  colour and marker of its value (with the markers' edge; a grey bar and a hollow circle without a value; a
  failed sample, listed without a bar, keeps its marker); without one, the colour and square of its group of
  identical genomes (a grey bar without a square in no group); without a comparison either, a single blue. The captions say so; the colour column's legend is the heatmap's, and a line under the bar charts
  gives it when there is no heatmap (no comparison, more than 150 genomes) or a value is held by failed samples
  only. With a colour column, the list of identical genomes gives each genome's marker and the samples table
  a hollow circle for a sample without a value; without one, the samples table starts the name of each genome
  of a group with the group's square. The tree's group squares hover as the heatmap's bands (`a: group 1`).
- The metadata column that colours the report may have up to 48 distinct values (was 8): each value gets a
  marker of its own, a colour and a shape, in the tree (instead of the circle), in the heatmap's band and legend
  (with its count), in the table of identical genomes by value and in the samples table. The 12 colours come
  from published colourblind-safe palettes (Okabe & Ito; Paul Tol's high-contrast, vibrant, muted, light and
  bright schemes), chosen and ordered so that every pair stays as distinct as possible under simulated
  deuteranopia, protanopia and tritanopia (Machado et al. 2009), the first ones, which a column with few values
  gets, the most; in the dark theme four colours change: the dark blue and the wine (under 3:1 on the dark
  surface) to a light cyan and a pink, the sky blue to Tol's bright yellow (not a second pale blue sharing the
  diamond with the cyan), and Tol's yellow to Okabe & Ito's (too close to the bright yellow). There are seven
  shapes (circle, triangle, square, diamond, inverted triangle, plus, cross). The first 12 values take the
  colours in turn with the first four shapes in turn, so that the values sharing a shape are four colours apart
  and stay distinct for a colourblind reader (at least 11 in CIEDE2000 and OKLab under the three simulated
  deficiencies, light and dark); the next 36 were chosen by a search to keep the values sharing a shape as
  distinct as possible for each number of values: as distinct as within the first 12 up to 25 values, then at
  least 9.4 up to 32, 8.6 up to 36, 7.6 up to 40 and 6.9 at 48 (taking four shapes in turn and moving them on
  after each round of 12 colours gave 5.9 from 13 values). So, up to 25 values, values whose colours look alike
  to a colourblind reader differ in shape; above, a few pairs sharing a shape may look alike, and the legend and
  the value written after each name tell them apart. Consecutive values differ in both colour and shape, and the
  markers follow the sorted values (adding or removing a value can change those of the values after it). A
  genome without a value keeps its hollow circle. With more than 12 values, and more values than groups, the
  table of identical genomes by value has the values as rows, so that it fits a printed page. The groups of
  identical genomes take their colours from the same palette, 12 instead of 8 before the grey, and keep their
  squares, now with the markers' thin dark edge (on the bands, the blocks and the list too: the first colour, a
  yellow, has 2:1 on the light surface); with a colour column they stay grey. A report rebuilt from a metadata
  copy edited by hand checks the recorded colour column again (more than 48 values: the figures are left
  uncoloured, with a note). The free-text rule is unchanged (from 10 samples with a value, a column with more
  distinct values than half of them is free text). As the limit was 8, a column with 9 to 48 values may now be
  the default colour column of a run made again (a report rebuilt from an existing run keeps its recorded one).
- The report's palette also colours the rest of the report: bars in blue (`#0077BB`), on the genome map the
  SNP ticks in blue where every genome has a call and in red where some genome has none, protein-coding genes
  in green, tRNA and rRNA genes in purple, and pseudogenes as a faint green fill with a green outline (3:1 on
  both surfaces).
- A gene without a symbol (no `/gene`, as in older RefSeq records: only `/locus_tag` and `/product`) is shown as
  `locus_tag (product)`, `LK299_pgr007 (23S ribosomal RNA)`, in the SNP table, the hovers, the intergenic
  contexts and the summary of the genes with the most SNPs (`LK299_pgr007 (23S ribosomal RNA; 4)`); on the map
  by its product when that is at most 12 characters (`tRNA-Val` for rice's `OrsajCt141`), else by its locus tag.
  Its identifier (the locus tag) is unchanged, and no symbol is made up from a product. In GFF3, a `Name` that
  is the locus tag, the product, the `gene_id` or the ID is not a symbol, nor the ID without its type prefix
  when the gene has a locus tag or the `Name` looks like one (`ID=gene-matK;Name=matK` keeps matK).
- GFF3: Ensembl's `ncRNA_gene` features are genes (their tRNA and rRNA transcripts are named after them, and a
  gene without a `Name` after its `gene_id`); a CDS or RNA feature without a `Parent` joins the gene feature
  containing it on its strand (the rRNA genes NCBI writes with `gene_biotype=other`: tomato NC_007898.3 now has
  the 142 genes of its GenBank record). GenBank: a quoted value continued on the next line is joined without a
  space after a hyphen, or after a comma followed by a digit, where the flat file broke a word
  (`2,6-diaminopimelate`).
- The genome map's caption says, for a reference with several sequences, that the band is above the genes or
  under the axis for the sequences without annotation, and the Methods paragraph describes the regions only when
  the report shows the band (not from `run_info.json` when no map was drawn).

### Fixed
- Paths with spaces (or other shell characters) anywhere (reads, reference, added genomes, output folder) made
  SKA2 (`ska build`), Parsnp, BBDuk and Flye fail. SKA2 now reads paths relative to its folder; Parsnp runs on
  links in `4_compared/.parsnp_work` (removed afterwards); BBDuk on links in `1_extracted`; Flye through links in
  a temporary folder (TMPDIR must not contain spaces).
- Moving a whole project (reads and output together) reran every step: an input file is now recognised by name,
  size and time when only its folder changed.
- Memory and CPU limits of cgroup v2 groups below the root (systemd, SLURM) and of cgroup v1 are detected; a CPU
  quota (rounded up) now limits the default `--threads`, with the CPU affinity.
- `--add-genomes` is checked before any step (a missing file says "file not found"; an invalid or used name is
  reported at once); with `--snp-method none` it is ignored with a warning.
- Metadata and sample sheets: UTF-16 and UTF-32 with a byte order mark are read; a file that is not UTF-8 is
  read as Windows-1252 with a warning, but a UTF-8 file with a few invalid bytes (a UTF-8 byte order mark, or
  any accented letter written in UTF-8) stays UTF-8, the invalid bytes read as `�` with a warning naming the
  first line (0.3.6 read the whole file as Windows-1252, so a path with `é` pointed to no file); only LF, CRLF and CR end a line. CSV is parsed as in 0.3.5 again (`"s1" ,file`),
  but an unclosed quote, also on the last line, is an error naming its line. Values in quotes keep their quotes in
  `metadata.tsv`; an existing `metadata.tsv` in another encoding no longer stops the run.
- An output path that is a file, or a folder that cannot be written, is a clear error; a checkpoint that cannot be
  saved (full disk) stops the run cleanly. Sample and added-genome names ending in a line break are refused. MD5
  works on FIPS hosts.
- SNP effects, GFF3: a spliced tRNA or rRNA written as one line from its first base to its last with `exon`
  children (NCBI: plastid *trnK*, *trnL*, *trnV*, *trnI*, *trnA*, *trnG*) has its exons only, so its intron is
  `intron`, as from the GenBank record (*matK*, in the intron of *trnK*, was `CDS / tRNA`); a part written
  wholly beyond the end of a circular sequence (the 5′ exon of *Epifagus* *rps12*, 82,777–82,890 on 70,028 bp)
  is moved back by the length instead of dropping the whole CDS with a "beyond the end" warning; a line across
  the origin among numbered parts no longer loses the other parts; a reverse-strand CDS listed by ascending
  coordinate is read downwards even when it spans more than half the sequence; and a trans-spliced CDS listed by
  ascending coordinate without `part=` numbers, in a sorted file (Ensembl), gives no effects, with a warning,
  instead of being read in a guessed order (potato's IRb *rps12* read exon 3 first). On seven NCBI records with
  both files (plastids, mitochondria), GenBank and GFF3 now give the same effects at every coding position and
  the same context at every position (but for a tRNA product the two files of yeast's mitochondrion name
  differently).
- SNP effects: a base read twice through a −1 ribosomal frameshift (`join(66..327,327..1228)`, F plasmid
  NC_002483.1; *E. coli* *dnaX*) changes every copy, one effect per codon (only the first copy was changed); a
  CDS with a part on another sequence (`join(X12345.1:1..100,201..500)`) gives no effects instead of effects in
  the wrong frame; a `/transl_except` codon split by an intron (`pos:join(...)`) is recognised; a sequence said
  to be a mitochondrion (`/organelle`, GFF3 `genome=`) without `/transl_table` uses table 1, the INSDC default
  (plant mitochondrial records give none), instead of 11; a translation table BACoN does not know (16, 21, 22,
  23...) gives a warning.
- Regions: a feature marking a junction or a border (`IRB/SSC junction`, `IRA-SSC border`), or any feature
  shorter than 50 bp, is no longer merged into an inverted repeat (IRb grew by 1 bp).
- Report: the Methods paragraph says that SNP effects were derived only when the genome map annotated SNPs (not
  without a VCF, a failed or skipped comparison, or an annotation that could not be used); the N-bases chart is
  replaced by a sentence for de novo runs (every bar was NA); the SNP table's Region column, and the summary's
  counts per region, include the SNPs of sequences without annotation, from their detected band; the N track
  sums the assemblies of the run's samples only (not those of removed or failed samples left in
  `3_assembled/all_assemblies/`); the genome map and the SNP table show the VCF records that passed their
  filters only (Parsnp writes the SNPs it left out of its alignment, and so of the distances, with `FILTER`
  `ALN`, `CID`, `LCB`...), and the caption counts the others; a printed figure fits one page (the heatmap was
  split); `MD5 ?` instead of `MD5 None`; the header says `report built with BACoN X` when that differs from the
  version that ran; `1 polishing iteration`; a deep tree (a
  caterpillar of 1,000 leaves) no longer exceeds Python's recursion limit in the report, `tree.svg` and the
  MultiQC files (every walk of a tree is iterative, and the MultiQC heatmap falls back to the table's order
  when the tree cannot be read); the MultiQC heatmap shows integers (`10`, not `10.00`) with recent MultiQC
  (1.35; 1.19 ignores the setting).
- Sample names may not start with a dot: a sample sheet naming a sample `..` deleted the output folder (and
  the reads in it), `.` every other sample's assemblies; the same rule applies to added genomes, and an assembler
  only deletes a sample's own folder in `3_assembled`.
- Parsnp 2.1.1 writes every SNP one base off (the reference never carries its own allele); BACoN dropped all of
  them and reported none. It now stops the comparison with an error, and requires Parsnp 2.1.2 or later.
- A comma in a path made Flye and IQ-TREE fail (both split their input on commas): Flye goes through its
  temporary folder, IQ-TREE runs in the comparison folder with relative paths.
- A BBDuk failure was always reported as "BBDuk ran out of memory" (the `java` line BBDuk prints names the
  out-of-memory option): only a real out-of-memory error says so. The log gives the working folder of a command
  run in one (`$ (cd <dir> && ...)`), so that it can be run again.
- A path BACoN may not read or create (permission denied) is a one-line error instead of a Python traceback. A
  linked `--add-genomes` file, or a single `-i` file, is named after the link, not its target. A failed or
  interrupted sample no longer leaves its merged reads (`1_extracted/.X.merged...`), minimap2's `.X.paf`, the
  `--keep-bam` read list or partial extracted reads behind. A checkpoint that cannot be saved no longer logs
  "Interrupted", and only a full disk or quota gets the "is the disk full?" hint.
- Annotation: NCBI GFF3 writes the inverted repeats of a plastome as `inverted_repeat` features, which were not
  recognised (the GFF3 of *Hydnora* NC_029358.1 gave no band, of *Lemna* NC_010109.1 one 2 bp off): GenBank and
  GFF3 now give the same regions on 18 NCBI records. GenBank: an RNA or CDS without a gene name only joins a gene
  containing it (an ncRNA overlapping the end of rrn26 made it an ncRNA gene, in the Arabidopsis mitochondrion
  NC_037304.1); a gene's kind is set by its CDS, then its tRNA or rRNA, then other RNAs, whatever their order; and
  a protein-coding gene's exons are its CDS parts (an ncRNA in the rps3 intron made 168 intron positions `CDS`).
  GFF3: in a sorted file without `part=`, a part written beyond the end of a circular sequence put the parts in
  the wrong order (wrong effects at 372 positions of rice's rps12); exon lines before their RNA's line were
  lost. A CDS with a stop codon inside its coding sequence (not from `/transl_except` or RNA editing) gives no
  effects, with a warning (a wrongly ordered CDS). `/organelle="mitochondrion:kinetoplast"` is a mitochondrion
  (table 1 by default). A mitochondrion's annotated inverted repeats no longer give a "not a plastome layout"
  warning.
- Trees: a deep tree of identical internal nodes (zero-length branches, as FastTree writes for identical
  genomes) no longer exceeds Python's recursion limit when rooted; a branch length or support that is not a
  finite number (`1e999`, `nan`) is an error; a tree that cannot be drawn loses its figure only, with a note,
  not the whole report.
- Report: the bar charts no longer shorten names of 4 or 8 characters (or 16 and 17 with markers) that fit;
  with more than 150 genomes the tree's caption points to the legend under the bar charts, which then also
  lists the reference and added genomes; a printed figure keeps its size (small figures filled the page); the
  header says that the last run's time does not count the steps it reused (`0.0 s` for a resumed run).
- Example: `check_example.py` prints the report check's OK line when that check passes, whatever the others.
- Docs: Methods and Outputs say that SNPs within about 30 bases of the reference's ends are not compared by
  Parsnp (15 by SKA2 with templated assemblies), what the `FILTER` values of a Parsnp VCF mean, and that the
  report shows the `PASS` records; the Usage table lists `-v/--version` and `-h/--help`; `CITATION.cff`'s
  top-level DOI is the concept DOI, so "Cite this repository" gives a DOI that resolves at any tag; the release
  steps say to rename the changelog section (a test checks it) and to take the wiki's report images again.

## 0.3.6 (2026-10-07)

### Added
- Reference annotations: `--annotation FILE` (GenBank or GFF3, gzipped or not), or a GenBank file as `-r`
  (its sequence is the reference, named as in NCBI's fasta of the record, and its features the annotation).
  The report's genome map then shows the genes of each strand, the LSC/IRb/SSC/IRa regions of a plastome
  derived from its inverted repeats, and labels for the genes with the most SNPs, and a new sortable table
  gives each SNP's region, gene, context and effect on the coding sequence (codon and amino-acid change,
  synonymous/missense/nonsense/stop lost or retained/start lost or retained; translation table 11 unless the
  annotation says otherwise, tables 1–5, 9, 11, 13 and 14 known; no effect at a `transl_except` codon). Genes
  in several pieces (trans-spliced, or across the origin of a circular sequence) are read from either format
  in their order of translation. The annotation is for the report only: it is not part of any checkpoint, so adding it to a
  finished run reruns nothing; it is copied to `OUTPUT/annotation.gb` or `.gff3` and recorded in
  `run_info.json`. `snps.vcf` is unchanged.
- Sample metadata: `--metadata FILE` (TSV or CSV with a `sample` column and any other columns), and the columns
  of a sample sheet other than `sample` and `file` (`--metadata` wins column by column). The report shows the
  columns in the samples table (sortable, numbers as numbers), and one column, `--color-by` or the first with at
  most 8 distinct values that is not free text (`none` for no colours), colours the tree (a circle and the
  value after each name) and the heatmap (a second band on both axes, with a legend), with a table of the groups
  of identical genomes against its values. The metadata columns are also in the MultiQC samples table. For the
  report only: not part of any checkpoint, copied to `OUTPUT/metadata.tsv` and recorded in `run_info.json`.

### Changed
- `report.html` has a new look (light and dark themes, numbered figures) and new figures: bar charts of the
  depth and of the `N` bases of each sample; the tree drawn with a coloured square for each group of identical
  genomes, the supports and a scale bar in substitutions per site with its equivalent in SNPs; a heatmap of the
  distances with colour classes fitted to their range and the groups of identical genomes as coloured bands;
  and a genome map with the positions of the VCF's SNPs and, for templated assemblies, the `N` bases along the
  reference.
- The bundled example is an annotated, plastid-like reference (LSC, IRb, SSC, IRa; 23 synthetic genes, as
  GenBank and fasta) with a sample metadata file; its 20 SNPs fall in chosen genes and contexts, and
  `example/run_example.sh` checks their effects as well as the distances, so the example shows every feature of
  the report. The published reports and the wiki pictures were rebuilt with all the features.
- The report's page is wider (1,440 px at most), so that a samples table with metadata columns fits; printed,
  its tables are complete (landscape, cells wrap).
- Sample sheets and metadata files are read strictly where it matters: an unclosed quote in a CSV, or a quoted
  value over several lines, is an error naming the line; a TSV's cells are split on tabs only (a cell entirely
  in quotes loses them, as before; other quotes are kept as written); two columns with the same name are an
  error saying so. In a sample sheet, lines starting with `#` are comments before and after the header, as
  before; in a metadata file, only before it (so a value such as `#FF0000` may come first).
- A moved or copied output folder resumes without rerunning anything, even when the copy did not keep the files'
  times (`cp -r` without `-a`, `scp`, an archive): in such a folder, the input of each step and of the
  comparison is recognized by its size (0.3.5 ran filtering, assembly and comparison again).
- The report's SNP table has a `Sequence` column when the reference has several sequences, and its summary
  counts the SNPs on sequences without annotation. The N track rescales a templated assembly to the reference
  only when their lengths differ by 5% or less (indels); a shorter or longer record is counted at its own
  positions. A `metadata.tsv` found in the output folder but not given to the run is said to be so in the
  report's metadata line (`run_info.json` and the MultiQC table do not have it).
- GFF3 annotations: the lines of a CDS across the origin listed by ascending coordinate (a sorted file) are read
  in the order of translation on both strands, with the sequence length from `##sequence-region` or, without
  it, from the reference; a single line ending beyond the length of a sequence flagged `Is_circular=true`
  (Bakta's way of writing a feature across the origin) is split into its two parts instead of being ignored.

### Fixed
- An unclosed quote in a sample sheet silently swallowed the rows after it (their samples were not run).
- A first run interrupted before its end left its copies `annotation.gb` (or `.gff3`) and `metadata.tsv`
  looking hand-made to the next run, which kept them when it had no annotation or metadata: the copies are now
  recorded in `.checkpoints/copies.json` as they are written, not only in `run_info.json` at the end.
- The annotation copy replaced non-ASCII characters with `?`; the file's bytes are copied as they are
  (decompressed when gzipped).
- `python -m bacon.report` on a copied output folder read the comparison files of the original folder through
  the absolute paths of `run_info.json`; it reads the folder's own files, and notes a missing one. On a run that
  did not finish, it says so (resume it to get a report) instead of "not a BACoN output folder".
- The report's VCF reader classed each distinct sample column instead of each distinct genotype, which was slow
  with per-sample fields (`GT:DP`...).
- The report: the VCF is read about 25 times faster for large runs; the support of a branch at the root could be
  drawn outside the tree; heatmap group labels overlapped with many small groups; long sample names were cut in
  the bar charts; the small grey text, the "no value" markers and the numbers in some heatmap cells had too
  little contrast.

### Upgrading to 0.3.6
- Sample sheets of 0.3.5 are read as before (quoted TSV cells, lines starting with `#` after the header), except
  that an unclosed quote in a CSV is an error instead of silently losing the rows after it.
- Output folders of 0.3.2 to 0.3.5 resume as they are; nothing runs again. A copy made without the files' times
  no longer reruns filtering, assembly and comparison.

## 0.3.5 (2026-10-02)

### Changed
- BACoN stopped by SIGTERM or SIGHUP exits with 143 or 129 (128 + the signal), Ctrl-C with 130.
- SNP alignments (`ska.snps.fasta`, `parsnp.snps.fasta`) keep only the columns with at least two different
  nucleotides: with `--ska-min-freq` below 1, columns made of one nucleotide and gaps were counted as SNP sites
  (more than 95% of them in a simulated test) and given to the tree program; with any setting, a column whose
  only variation was an ambiguous base in one genome was counted too. The SNP distances do not change; the
  number of SNP sites and the trees' input do. Comparisons are not redone on resume: use `--redo compare`.
- The default `-t` and `-m` are the CPUs and memory BACoN may use under a job scheduler, `taskset` or in a
  container, not those of the whole machine.
- A templated assembly notes the reference sequences that no read covers (they are not in the assembly).

### Fixed
- VCF from Parsnp: positions where the reference's base is N (an ambiguity code in the reference) were written as
  a SNP of every genome.
- Midpoint rooting dropped the support of the split on which the root is placed (often the deepest one) from
  `tree.nwk`, `tree.svg` and the report.
- An output folder inside the input folder (`-i reads -o reads/bacon`) worked once, then every resume failed with
  an unrelated message (its files were taken for a sample); it is now refused.
- `-b bbduk`: a resumed run gave each sample its share of `-m` as if all samples were baited at the same time.
- A sample with no read baited lost its numbers of raw reads and bases in `summary.tsv`.
- Files in hidden folders inside a sample's folder were used.
- Read names with non-ASCII characters could not be matched, so their reads were not baited.
- A checkpoint of 0.3.3 or earlier with a sample named `distances` stopped the run.
- The report: the reference's file name was not escaped; the tree of a Parsnp comparison was said to be built on
  SNPs (it is built on the core-genome alignment); browsers that darken pages could make the tree invisible.
  MultiQC: the heatmap was said to be in tree order when there was no tree.
- `python -m bacon.report` on a folder that is not a BACoN output gave a traceback.
- `run()` left the `bacon` logger at the INFO level.
- `nohup bacon ... &` stopped when the terminal was closed (0.3.4 replaced nohup's ignored SIGHUP with its own
  handler); an ignored SIGHUP or SIGTERM now stays ignored.
- A resumed output folder of 0.3.3 or earlier could still reuse a result made from an earlier input (new reads
  baited, then a run interrupted before filtering): such results are now reused only if newer than their input,
  and get the input signatures of 0.3.4.
- The final error of a run (all samples failed, a failed comparison, "Interrupted") was not written to
  `bacon.log`.
- `bacon.cli.main` failed when called from a thread other than the main one, and left its signal handlers
  installed after returning.
- An interruption waited for the reads being copied or counted by BACoN itself (not by a program) to finish.
- A moved output folder of 0.3.3 or earlier whose path contained the name of a BACoN folder (such as
  `/data/1_extracted/out`) ran every step again.
- The upgrade notes: with `-b bbduk`, every step runs again once after an upgrade from 0.3.2, not only baiting.

## 0.3.4 (2026-10-01)

### Changed
- A VCF written by BACoN 0.3.2 or earlier is rewritten on resume, from the comparison's files (no need for
  `--redo compare`).
- SIGTERM and SIGHUP (`kill`, a closed terminal, a scheduler's time limit) stop BACoN like Ctrl-C: the programs
  it started are stopped. In 0.3.3 they kept running.

### Fixed
- 0.3.3 ran every step again on output folders of 0.3.2 (its checkpoints differed even with minimap2).
- A sample's result could be reused after its input changed, when a run was interrupted in between: new reads
  baited, then the run stopped before filtering; or a step run with other parameters, interrupted, and run
  again with the first ones. Each result now records the input it was made from, and so does the comparison.
- A moved or copied output folder resumed with the files of the original folder.
- After an interruption, the programs killed were reported as failures of their samples.
- On a file system without locks (some NFS or SMB mounts), BACoN said another run was using the folder.
- The report's methods said BBDuk allowed two mismatches, whatever `--hdist`.

### Upgrading to 0.3.4
- From 0.3.2: output folders resume, and their VCF is rewritten. Steps run again once only for a soft-masked or
  IUPAC reference (every step), `--add-genomes` (the comparison) and `-b bbduk` (every step, as baiting runs
  again for the new `--hdist` default).
- From 0.3.3: output folders baited with minimap2 run every step again once (0.3.3 changed the checkpoints by
  mistake); with `-b bbduk`, they resume.

## 0.3.3 (2026-10-01)

### Added
- `--hdist`: the mismatches BBDuk allows in a k-mer; the default is now 1 (was 2), which baits the same reads
  in the example with a fraction of the memory. BBDuk running out of memory gives a hint.

### Fixed
- An interrupted run (Ctrl-C) lost the samples finished in the step it was in; each sample is now recorded
  when it finishes. The programs still running, and their children (Flye's minimap2), are stopped.
- Two runs in the same output folder at the same time overwrote each other's files; the second one now stops.
- The traceback of an unexpected error in a sample is written to that sample's log.
- A file listed twice for one sample in a sample sheet doubled its reads; it is used once, with a warning.
- An `--add-genomes` genome could take the name of a sample that failed.
- SKA2: added genomes counted towards the rule that decides whether the genomes are circular, and were not
  extended when the reference was; the rule now counts the samples' assemblies, and the added genomes are
  extended with the reference.
- `--keep-bam` read all the read names through Python; `summary.tsv`, the distances, the VCF and the MultiQC
  files are written atomically; the report shows no inline heatmap above 150 genomes (the file is linked).
- A soft-masked (lower-case) reference lost the SNPs of its lower-case regions from `snps.vcf` (`ska map` and
  HarvestTools are case-sensitive), and IUPAC codes in a reference or added genome were read by SKA2 as fixed
  bases (false SNPs): the reference, the assemblies and the added genomes are now upper case with N for any
  ambiguity code.
- An unexpected error in one sample (a truncated or corrupt gzip file, a non-ASCII read name) stopped the whole
  run without a summary; it now fails that sample only, with a clear message. An unexpected error in the
  comparison still writes the summary, `run_info.json` and the report.
- `--add-genomes` made the comparison run again on every resume.
- Reads from fasta files were all kept when none reached `--min-read-length`.
- A failure in a pipe of programs was attributed to the upstream program killed by SIGPIPE.
- Parsnp: a warning when the core genome is less than half of the reference (a short genome can shrink it for
  every genome).
- Switching from a de novo to the templated assembly left the de novo assembly graphs.
- The VCF of a comparison made by BACoN 0.3.1 or 0.3.2 now uses the reference as SKA2 used it (extended when circular).
- Genotypes such as `0/1` or `1|0:35` are renumbered correctly when an `N` allele is removed.
- `example/run_example.sh` checks the SKA2 comparison, and its output folder is ignored by git.
- VCF: `N` (ambiguous or missing bases, such as the N of templated assemblies) was written as an alternate
  allele, and Parsnp wrote one record per `N`; it is now a missing genotype, and records left without an
  alternate allele are dropped, as are records where the reference's own base is ambiguous.
- VCF from SKA2: SNPs within 15 bp of the ends of circular genomes were in the alignment but not in the VCF;
  `ska map` now uses the reference as SKA2 used it (`ska_reference.fasta`), and positions are folded back.
- VCF from Parsnp: HarvestTools' undeclared INFO value `NA` made bcftools warn, and `bcftools norm` fail.
- A comparison made by BACoN 0.3.1, or whose VCF failed, never got a VCF on later runs; it now does, from the
  comparison's files, without redoing the comparison.
- A failed `ska map` left `snps.raw.vcf`; a malformed VCF from the tools stopped the run or gave an invalid
  file: it is now reported and the run goes on.

### Upgrading from 0.3.2
- Every step runs again once: the checkpoints changed (by mistake; 0.3.4 resumes folders of 0.3.2, with the
  exceptions below).
- A soft-masked or IUPAC reference rewrites `reference.fasta`, so every step runs again once.
- Comparisons with `--add-genomes` are redone once.
- `-b bbduk`: every step runs again once, as baiting runs again for the new `--hdist` default.

## 0.3.2 (2026-09-26)

### Added
- `snps.vcf` in the comparison folder: the SNPs of every genome relative to the reference (SKA2: `ska map`;
  Parsnp: HarvestTools), with contig lengths in the header, one column per genome
  ([#1](https://github.com/duceppemo/BACoN/issues/1)).

### Fixed
- MultiQC: the files of several runs in one search path were merged (samples pooled, only one heatmap kept);
  the sections are now named after the output folder.
- IQ-TREE renamed samples with `+` in their name, which then did not match in the report and MultiQC heatmaps;
  the tree programs now see placeholder names and the tree gets the sample names back.
- The report counted genomes as identical when they were connected through zero distances only (distances skip
  positions with `N`), and claimed identity when no SNP site was compared.
- The reference MD5 in `run_info.json` and the report was that of BACoN's copy, not of the file given.
- The methods paragraph of the report: fasta reads (no Filtlong, no quality cut), `--template-gaps
  reference`, Flye's genome size and overlap, added genomes, failed comparisons, IQ-TREE's consensus tree,
  Parsnp's core-genome tree, and the circular-contig extension only when there are circular contigs.
- `python -m bacon.report` on a moved output folder lost the tree and distances; incomplete `run_info.json`
  or an empty `summary.tsv` no longer crash it.
- Tabs or line breaks in an error message could shift the columns of `summary.tsv`.
- The report keeps its colours when printed, and sorts text columns naturally.

## 0.3.1 (2026-09-26)

### Added
- `report.html`: a self-contained HTML report written at the end of every run: overview, sortable sample table
  with flagged values, tree, SNP-distance heatmap in tree order, identical genomes, a methods paragraph with
  program versions, and provenance. `python -m bacon.report OUTPUT` rebuilds it.
- MultiQC custom content: `bacon_samples_mqc.json` (table), `bacon_reads_mqc.json` (bases kept, filtered out
  and off-target), `bacon_distances_mqc.json` (SNP-distance heatmap).
- `run_info.json` records the reference file, number of sequences, length and MD5.

## 0.3.0 (2026-09-25)

A rewrite. Every program was reconsidered: those no longer maintained were replaced after comparing the
candidates on simulated data of known truth and on real data ([Validation](https://github.com/duceppemo/BACoN/wiki/Validation)).

### Changed
- Installable package with a `bacon` command (`pip install .`); `python bacon.py` still works. BACoN itself
  uses the Python standard library only (pandas, psutil, pysam and ete3 are no longer needed).
- The default assembly is now **templated**: the consensus of the reads aligned to the reference, with
  `samtools consensus` (`-a samtools`). It replaces Rebaler and was the most accurate method in every test.
  Flye remains available (`-a flye`), and myloasm is added (`-a myloasm`).
- The default SNP method is now **SKA2** (`--snp-method ska`), reference-free, exact on all simulated data and
  counting a SNP in an inverted repeat once. `--ska-min-freq` below 1 gives pan-genome SNPs, the role of kSNP3
  in BACoN 0.1. Parsnp remains (`--snp-method parsnp`), now with `-c`, which turns off its filter that drops divergent genomes.
- Output folders: `1_extracted`, `2_filtered`, `3_assembled`, `4_compared/<method>/` (there is no trimming
  step any more). Each comparison method writes to its own folder.
- Resuming is automatic and parameter-aware: finished steps are skipped when their parameters and inputs did
  not change; changing a parameter reruns that step and the following ones; failed samples are retried;
  adding samples processes only them (and the comparison). `--redo STEP` forces a step. The `done_*` files
  are replaced by `.checkpoints/`.
- Trees are midpoint-rooted and drawn as SVG by BACoN (`tree.svg`), without Qt.
- Sample names are the file names without their extension (`_pass` is no longer removed); a folder of files
  (such as MinKNOW's `barcode01/`) is one sample.
- `--kmer-size` defaults to 31, the largest value BBDuk accepts (99 made BBDuk fail).

### Added
- `summary.tsv`: reads and bases at each step, baited share, depth, assembly statistics, circular contigs,
  status and notes (low depth, unexpected length, `N` bases, failures) for every sample.
- `snp_distances.tsv`: pairwise SNP distances.
- `run_info.json`: command, settings, versions of every program, samples and results; `bacon.log`; per-sample
  logs of every program in `logs/`.
- `--add-genomes`: finished genomes (such as published plastomes) included in the comparison.
- `--sample-sheet` (several files per sample), `--tree iqtree`, `--read-type`, `--template-gaps`,
  `--target-depth`, `--keep-percent`, `--min-read-length`, `--debug`.
- A failed sample no longer stops the run; every program's exit status is checked.
- Tests (pytest, stub programs), continuous integration, a simulated example with known SNPs
  (`example/`), a validation suite (`validation/`), a bioconda recipe, and a wiki maintained in `docs/wiki/`.

### Removed
- Porechop (unmaintained; reads are trimmed by the basecaller), Shasta, Rebaler, Snippy (cannot be installed
  with current assemblers), PhaME (its version check rejects samtools 1.10 to 1.29), RAxML and the
  `requirements.txt` environment. See [Methods](https://github.com/duceppemo/BACoN/wiki/Methods#choices-that-changed-in-03).

### Fixed
- With minimap2 baiting, the BAM files of other samples processed in parallel could be deleted before their
  reads were extracted.
- `--threads` and `--parallel` were not checked (the bounds test could never be true).
- Snippy ran with a changed working directory shared by parallel threads.
- PhaME added the reference to the assembly folder, which later runs took for a sample.
- Shasta needed `pigz`, which was not in the environment.
- Programs failed silently: their output was discarded and exit codes ignored.
- Reads given as fasta failed at filtering (Filtlong needs qualities); they are now filtered by length.

## 0.2 (2023-03-20)

- Snippy (default) and PhaME as SNP methods; kSNP3 removed.

## 0.1 (2022-08-30)

- First release: baiting (minimap2, BBDuk), Porechop, Filtlong, Flye/Shasta/Rebaler, Parsnp and kSNP3.
