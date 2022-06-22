# NanOrgAssM
NanOrgAssM stands for "Nanopore Organelle Assembly Method". As its name suggests, it assembles organelle genomes from Nanopore (genome skimming) data. Taking advantage of targeted sequencing to enrich for organelle sequences is recommended, especially for larger genomes (I like to call that "Power Skimming").

## Workflow
NanOrgAssM follows these steps:
1. Extract organelle Nanopore fastq reads with `minimap2` using a provided reference sequence (ideally the same used for the targeted sequencing).
2. Remove Nanopore adapters with `porechop` and filter low quality extracted reads with `filtlong`.
3. Assemble the extracted reads with `flye`. Note that the `ont-hq` flag is used here, assuming you used Guppy5+ or Q20 (<5% error).
4. Compare samples with a core SNP approach usinf `parsnp`

## Installation
1. Create virtual environment and install all dependencies. Requires conda to be installed. See [here](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links) for instructions if needed:
```
# Create environment
conda create -n nanorgassm -y -c bioconda -c etetoolkit porechop filtlong minimap2 \
    samtools flye bbmap git ete3 pysam bandage parsnp

# Activate enviroment
conda activate nanorgassm
```
2. Clone and test NanOrgAssM
```
# Clone repo
git clone 
```

## Example
1. First step is to basecall the fast5 using guppy using the terminal in Linux. You'll need to adjust the parameters based on your graphics card power, the accuracy wanted and if samples were multiplexed or not.
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
2. Run NanOrgAssM to extract and assemble organelle reads matching your reference.
```
# Change path according to where you cloned NanOrgAssM on your machine.
# Don't forget to activate your nanorgassm virtual environment if not already done.
python nanorgassm.py \
  -r /path/to/my_reference.fasta \
  -i /path/to/folder/containing/Nanopore/reads \
  -o /path/to/nanorgassm/output/folder \
  -t $(nproc)
  
```
`-r` is the reference sequence used to bait the reads
`-i` is the input reads. Fastq and fasta are supported, gzipped or not.
`-o` is the output folder.
`-t` is the number of threads used. `$(nproc)` means to use all availble threads (and should not be used on shared computing resources).
usage: nanorgassm.py [-h] -r /path/to/reference_organelle/genome.fasta -i /path/to/input/folder/ -o /path/to/output/folder/ [-t 16] [--parallel 2]

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
  --parallel 2, -p 2    Number of samples to process in parallel. Default is 2. Optional

## Usage
```commandline
usage: python nanorgassm.py [-h] -r /path/to/reference_organelle/genome.fasta -i /path/to/input/folder/ -o /path/to/output/folder/ [-t 16] [-p 2]

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

## Cheats
As it's finishing the various steps, check point files are created so NanOrgAssM can resume where it left off is something wrong happened during the run. In order:
```commandline
done_extracting
done_trimming
done_filtering
done_assembling
done_comparing
```
This can be exploited to combine multiple NanOrgAssM outputs and run only the last part to get the tree. The steps would be to:
1. Create a new input folder with all the assemblies to include (created with the reference sequence preferably).
2. Create a new output folder and create inside that folder the `done_extracting`, `done_trimming`, `done_filtering` and `done_assembling` files (e.g. `touch done_extracting`).
3. Run NanOrgAssM using the newly populated input and output folders.
