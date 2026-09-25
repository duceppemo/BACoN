# Tutorial: plastomes of 28 potato cultivars

This tutorial runs BACoN on public data: whole-genome Nanopore reads of 28 potato (*Solanum tuberosum*)
cultivars from the Ural region, sequenced to compare their plastomes (NCBI BioProject
[PRJNA807056](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA807056); MinION, reads longer than 6 kb, 100–700 Mb
per sample). Potato cultivars carry a few types of plastome ("cytoplasm types"); the common T-type differs from
the others by, among other things, a 241 bp deletion between the *ndhC* and *trnV-UAC* genes. Let's see what
BACoN finds.

## 1. Download the reads and the reference

About 7 GB. The script downloads the 28 runs from ENA, named after their cultivar, and checks their MD5:

```bash
mkdir potato && cd potato
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA807056&result=read_run&fields=sample_alias,fastq_ftp,fastq_md5&format=tsv" \
    | tail -n +2 > runs.tsv
mkdir -p reads
while IFS=$'\t' read -r run alias url md5; do
    name=$(echo "$alias" | tr ' ' '_')
    curl -s -C - --retry 5 -o "reads/$name.fastq.gz" "https://$url"
    echo "$md5  reads/$name.fastq.gz"
done < runs.tsv > md5.txt
md5sum -c md5.txt

# The reference: the plastome of cultivar Désirée (T-type)
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_008096.2&rettype=fasta" \
    > NC_008096.2.fasta
```

If an MD5 check fails, run the loop again: `curl -C -` resumes interrupted downloads.

## 2. Run BACoN

```bash
bacon -r NC_008096.2.fasta -i reads/ -o bacon_potato -t 32 -p 8
```

The defaults: baiting with minimap2, Filtlong capping each sample at 100x, templated assembly with samtools,
core SNPs with SKA2, and a FastTree tree. It takes about 3 minutes.

## 3. Reads and assemblies: `summary.tsv`

```
Sample     Status  Raw_reads  Raw_bases  Baited_reads  Baited_pct  Filtered_reads  Est_depth  Assembly_length  N_bases
12_22_134  ok      18388      202618009  3553          17.437      1433            100.0      155386           133
14_4_1     ok      15023      182354649  3018          18.349      1313            100.1      155383           153
14_6_3     ok      20396      234625686  4122          17.970      1292            100.2      155166           0
...
```

- 10–23% of the bases are plastid reads, typical of leaf DNA; three runs (Alaska, Argo, Shah) were
  already filtered to plastid reads (100%).
- Every sample reaches 65–100x after filtering, far more than needed.
- The consensus lengths fall into two groups: about 155,170 bp and about 155,390 bp. Samples of the second group
  have 120–150 `N` bases, mostly in one region, next to position 52,580.

## 4. SNPs and tree: `4_compared/ska/`

`snp_distances.tsv` gives 119 SNP sites between the 28 plastomes and the reference. The cultivars fall into
three groups:

| Group | Cultivars | SNPs to the reference | Within the group |
|---|---|---|---|
| T-type | 14_6_3, 16-35-5, 16_1_2, Bagira, Bankir, Iskra, Luks, Shah, Terra, Zdraven | 1 | identical |
| Lineage A | 15-27-1, Legenda | 67 | identical |
| Lineage B | 12_22_134, 14_4_1, 15_22_4, 16_4_3, Alaska, Amur, Argo, Baron, Bravo, Gornyak, Irbitskiy, Kamenskiy, Mishka, Otrada, Start, Utro_ranneye | 67–76 | 0–12 |

Lineages A and B are 66 and 75 SNPs from the T-type group, and 84 from each other.

![Tree of the 28 potato plastomes](https://raw.githubusercontent.com/duceppemo/BACoN/main/docs/images/tutorial_potato_tree.svg)

The T-type group shares the reference's plastome but for one SNP; lineages A and B carry other cytoplasm
types.

## 5. The 241 bp marker

The `N` bases of lineages A and B sit in the *ndhC*–*trnV-UAC* spacer (reference positions 51,834 to
52,696). Aligning an assembly of each group to the reference shows why:

```bash
minimap2 -cx asm5 --cs NC_008096.2.fasta bacon_potato/3_assembled/all_assemblies/Alaska.fasta
```

the plastomes of lineages A and B have an insertion of about 240 bp at position 52,578 relative to the T-type reference:
the known T-type deletion, seen from the other side. The templated assembly places it, but the reads do not
agree on all of its bases, hence the `N`. A de novo assembly resolves it completely:

```bash
bacon -r NC_008096.2.fasta -i reads/ -o bacon_potato -a flye -t 48 -p 12
```

reuses the baited and filtered reads and assembles each sample with Flye (about 10 minutes). The Flye
assembly of Alaska has an insertion of exactly 241 bp at position 52,578. Flye also reports 17 of the 28
plastomes as one circular contig. Its SNP distances (`4_compared/ska/snp_distances.tsv`, replaced) are
identical to those of the templated assembly for all 406 pairs of samples.

## What this shows

- Genome skimming data from a few hundred megabases per sample is enough for complete plastomes.
- The templated assembly gives the SNPs and small indels in minutes; a de novo assembly resolves the insertions
  and the structure, and confirms the SNPs independently.
- Identical plastomes (distance 0) are common among cultivars: the plastome is inherited maternally, and
  cultivars descend from few maternal lineages.
