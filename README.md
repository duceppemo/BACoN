# BACoN
BACoN stands for "Bait, Assemble and Compare Nanopore". It's extracts and assembles Nanopore reads matching a reference from many samples in parallel. Assemblies are then compared by extracting core SNPs and building a tree.

BACoN was designed with genome skimming in mind. BACon alows to reconstitute organelle genomes from Nanopore data.  Taking advantage of targeted sequencing to enrich for organelle sequences is recommended, especially for larger genomes (I like to call that "Power Skimming").

## TODOs
* Add a few `try` to handle potential errors: No reads baited. No assembly created (not enough reads or reads too short).
* Maybe provide a library of organelle genomes or auto-detect species with Mash from pre-compile organelle library and auto-download reference if not provided.

## Workflow
BACoN follows these steps:
1. Extract Nanopore reads matching reference with `minimap2` (ideally using the same reference used for targeted sequencing).
2. Remove Nanopore adapters with `porechop` and filter low quality extracted reads with `filtlong`.
3. Assemble the extracted reads with `flye`. Note that the `ont-hq` flag is used here, assuming you used Guppy5+ or Q20 (<5% error).
4. Compare samples with a core SNP approach using `parsnp`

## Installation
1. Create virtual environment and install all dependencies. Requires conda to be installed. See [here](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links) for instructions if needed:
```
# Create environment
conda create -n BACoN -y -c bioconda -c etetoolkit porechop filtlong minimap2 \
    samtools flye bbmap git ete3 pysam bandage parsnp harvesttools raxml fasttree

# Activate enviroment
conda activate BACoN
```
2. Clone and test BACoN
```
# Clone repo
git clone https://github.com/duceppemo/BACoN

# Test BACoN
cd BACoN
python bacon.py -h
```

## Usage
```commandline
usage: python bacon.py [-h] -r /path/to/reference_organelle/genome.fasta -i /path/to/input/folder/ -o /path/to/output/folder/ [-t 16] [-p 2]

Organelle genome assembly from Nanopore genome skimming data.

optional arguments:
  -h, --help            show this help message and exit
  -r /path/to/reference_organelle/genome.fasta, --reference /path/to/reference_organelle/genome.fasta
                        Reference genome for read mapping. Mandatory.
  -i /path/to/input/folder/, --input /path/to/input/folder/
                        Folder that contains the read files. Mandatory.
  -o /path/to/output/folder/, --output /path/to/output/folder/
                        Folder to hold the result files. Mandatory.
  -t 16, --threads 16   Number of threads. Default is maximum available(16). Optional
  -p 2, --parallel 2    Number of samples to process in parallel. Default is 2. Optional
```

## Example
1. First step is to base-call the fast5 using guppy using the terminal in Linux. You'll need to adjust the parameters based on your graphics card power, the accuracy wanted and if samples were multiplexed or not.
```
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
2. Run BACoN to extract and assemble organelle reads matching your reference.
```
# Change path according to where you cloned BACoN on your machine.
# Don't forget to activate your BACoN virtual environment if not already done.
python bacon.py \
  -r /path/to/my_reference.fasta \
  -i /path/to/folder/containing/Nanopore/reads \
  -o /path/to/BACoN/output/folder \
  -t $(nproc)
```
`-r` is the reference sequence used to bait the reads
`-i` is the input reads. Fastq and fasta are supported, gzipped or not.
`-o` is the output folder.
`-t` is the number of threads used. `$(nproc)` means to use all availble threads (and should not be used on shared computing resources).

3. Results
Here's a tree representation of the main output files and folders:
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
* `1_extracted` contains the baited reads by `minimap2`.
* `2_trimmed` contains the trimmed reads by `porechop`.
* `3_filtered` contains the filtered reads by `filtlong`.
* `4_assembled` contains the assembly sub-folders for each sample. All assemblies are located in the `all_assemblies` folder. Assembly quality can be assessed quickly from the `assembly_graphs` sub-folder.
* `5_compared` contains the SNPs and tree files from `parsnp`, `RAxML` and `FastTree`. RAxML and FastTree are provided to add bootstrap values to the tree.
* The `done` files are checkpoint files. See the `Cheat` section below for details.

## Cheats
As it's finishing the various steps, check point files are created so BACoN can resume where it left off is something wrong happened during the run. In order:
```commandline
done_extracting
done_trimming
done_filtering
done_assembling
done_comparing
```
This can be exploited to combine multiple BACoN outputs and run only the last part to get the tree. The steps would be to:
1. Create a new input folder with all the assemblies to include (created with the reference sequence preferably).
2. Create a new output folder and create inside that folder the `done_extracting`, `done_trimming`, `done_filtering` and `done_assembling` files (e.g. `touch done_extracting`). Note that we skipped creating the `done_comparing` file.
3. Run BACoN using the newly populated input and output folders.
