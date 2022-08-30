# BACoN
BACoN stands for "<ins>B</ins>ait, <ins>A</ins>ssemble and <ins>Co</ins>mpare <ins>N</ins>anopore". It extracts and assembles Nanopore reads matching a reference from many samples in parallel. Assemblies are then compared by extracting SNPs and building a tree.

BACoN was designed with genome skimming in mind. BACoN allows reconstruction of organelle genomes from Nanopore data.  Performing targeted sequencing to enrich sequence(s) of interest (say chloroplast genome) is highly recommended to achieve higher coverage, thus better assemblies, especially for larger genomes (I like to call this "Power Skimming").

## Workflow
BACoN follows these steps:
1. Extract Nanopore reads matching reference with `minimap2` or `bbduk` (ideally using the same reference used for targeted sequencing).
2. Remove Nanopore adapters with `porechop` and filter small (min 500bp) and lower quality (discard bottom 5% or keep top 100X) reads with `filtlong`.
3. Assemble reads with `flye`, `shasta` or `rebaler`. Note that the "ont-hq" flag is used for Flye, assuming that Guppy5+ or Q20+ chemistry (<5% error) was used for basecalling. Flye and Shasta are *de novo* assemblers while rebaler is a reference-guided assembler. Flye seems to work better with circular genomes. Shasta works well too with circular references (like full organelle genomes), but will also work with linear references (like ribosomal DNA). Rebaler is more experimental (not fully validated for this type of application, but I thought it would be good to have a "different" type of assembler to try in case the *de novo* methods don't work well).
4. Compare samples (make tree) with `parsnp` (core) or `kSNP3` (pan and core). Boostrap trees are also created using `RAxML` and `FastTree` from SNPs from parsnp.

## Installation
Requires conda. See [here](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links) for instructions if needed.
1. Clone repository:
```commandline
# Make sure you have "git" install in your base enviroment:
conda install git

# Clone repo to your prefered location:
git clone https://github.com/duceppemo/BACoN
```
2. Create virtual environment and test BACoN:
```commandline
# Create BACoN virtual environment:
cd BACoN
conda create -n BACoN -y -c bioconda -c etetoolkit -c hcc --file requirements.txt

# Activate enviroment:
conda activate BACoN

# Test program. No errors should be displayed:
python bacon.py -h
```

## Usage
```commandline
usage: python bacon.py [-h] -r /path/to/reference_organelle/genome.fasta -i /path/to/input/folder/ or /path/to/my_fastq.gz -o /path/to/output/folder/ [-b {minimap2,bbduk}] [-a {flye,shasta,rebaler}] [--min-size 3000] [-t 16]
                [-p 2] [-m 57] [-s 150000] [-k 99] [--keep-bam] [-snp {parsnp,ksnp3}]

Extract, assemble and compare Nanopore reads matching a reference sequence.

optional arguments:
  -h, --help            show this help message and exit
  -r /path/to/reference_organelle/genome.fasta, --reference /path/to/reference_organelle/genome.fasta
                        Reference genome for read mapping. Mandatory.
  -i /path/to/input/folder/ or /path/to/my_fastq.gz, --input /path/to/input/folder/ or /path/to/my_fastq.gz
                        Folder that contains the fastq files or individual fastq file. Mandatory.
  -o /path/to/output/folder/, --output /path/to/output/folder/
                        Folder to hold the result files. Mandatory.
  -b {minimap2,bbduk}, --baiting-method {minimap2,bbduk}
                        Read baiting method. Default "minimap2". Optional.
  -a {flye,shasta,rebaler}, --assembly-method {flye,shasta,rebaler}
                        Assembly method. Default "flye". Optional.
  --min-size 3000       Minimum read size for Shasta assembler or minimum read overlap for Flye. Default 3000 for Shasta and auto for Flye. Optional.
  -t 16, --threads 16   Number of threads. Default is maximum available(16). Optional.
  -p 2, --parallel 2    Number of samples to process in parallel. Default is 2. Optional.
  -m 57, --memory 57    Memory in GB. Default is 85% of total memory (57)
  -s 150000, --size 150000
                        Override automatically detected reference size. Optional.
  -k 99, --kmer-size 99
                        Kmer size for baiting. Only used if "--baiting-method" is "bbduk". Optional.
  --keep-bam            Do not delete BAM files. Only used if "--baiting-method" is "minimap2". Optional.
  -snp {parsnp,ksnp3}, --snp-method {parsnp,ksnp3}
                        SNP calling method. Default "parsnp". Optional.
```

## Example
1. First step is to perform basecalline (fast5->fastq) using the lastest `guppy`. You'll need to adjust the parameters based on your graphics card, the accuracy wanted and if samples were multiplexed or not.
```commandline
fast5=/path/to/folder/holding/fast5
basecalled=/path/to/folder/holding/fastq

# Basecall using RTX A4000 and Titan RTX.
# Q20+ kit using barcodes
guppy_basecaller \
    --config dna_r10.4_e8.1_sup.cfg \
    --input_path "$fast5" \
    --save_path "$basecalled" \
    --calib_detect \
    --recursive \
    --records_per_fastq 0 \
    --compress_fastq \
    --disable_pings \
    --gpu_runners_per_device 2 \
    --chunk_size 1000 \
    --chunks_per_runner 128 \
    --device "cuda:0 cuda:1" \
    --barcode_kits "SQK-NBD112-24" \
    --detect_barcodes \
    --detect_adapter \
    --detect_primer \
    --trim_barcodes \
    --trim_adapters \
    --trim_primers
```
2. Run BACoN to extract and assemble reads matching your reference.
```commandline
# Change path according to where you cloned BACoN on your machine.
# Don't forget to activate your BACoN virtual environment if not already done.
python bacon.py \
  -r /path/to/my_reference.fasta \
  -i /path/to/folder/containing/Nanopore/reads \
  -o /path/to/BACoN/output/folder
```
- `-r` is the reference sequence used to bait the reads.
- `-i` is the input reads. Fastq and fasta are supported, gzipped or not.
- `-o` is the output folder.

3. Results
Here's a tree representation of the main output files and folders using default settings.
```commandline
bacon
├── 1_extracted
│        ├── sample1.fastq.gz
│        ├── sample2.fastq.gz
│        └── sample3.fastq.gz
├── 2_trimmed
│        ├── sample1.fastq.gz
│        ├── sample2.fastq.gz
│        └── sample3.fastq.gz
├── 3_filtered
│        ├── sample1.fastq.gz
│        ├── sample2.fastq.gz
│        └── sample3.fastq.gz
├── 4_assembled
│        ├── all_assemblies
│        ├── assembly_graphs
│        ├── sample1
│        ├── sample2
│        └── sample3
├── 5_compared
│        ├── fasttree.pdf
│        ├── fasttree.tree
│        ├── parsnp.fasta
│        ├── parsnp.tree
│        ├── parsnp.xmfa
│        ├── RAxML_bestTree.raxml.tree
│        └── RAxML_bipartitionsBranchLabels.raxml.tree
├── done_assembling
├── done_comparing
├── done_extracting
├── done_filtering
└── done_trimming
```
- `1_extracted` contains the baited reads by `minimap2`.
- `2_trimmed` contains the trimmed reads by `porechop`.
- `3_filtered` contains the filtered reads by `filtlong`.
- `4_assembled` contains the assembly sub-folders for each sample. All assemblies are located in the `all_assemblies` folder. Assembly quality can be assessed quickly from the `assembly_graphs` sub-folder.
- `5_compared` contains the SNPs and tree files from `parsnp`, `RAxML` and `FastTree`. RAxML and FastTree are provided to add bootstrap values to the tree.
- The `done` files are checkpoint files.
