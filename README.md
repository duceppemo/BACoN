# NanOrgAssM
NanOrgAssM stands for "Nanopore Organelle Assembly Method". As its name suggests, it assembles organelle genomes from Nanopore genome skimming data. Hopefully, targeted sequencing with a suitable organelle reference sequencing has been used to enrich organelle sequences in the reads (I like to call that "powerskimming").

## TODO

## Information
NanoOrgAssM follows these steps:
1. Basecall and demultiplex the fast5 with `guppy` in SUP mode to generate the fastq files.
2. Remove adapters with `porechop` and filter low quality reads with `filtlong`.
3. Extract organelle reads with `minimap2` using a provided reference sequence (ideally the same used for the targeted sequencing).
4. Assemble the extracted reads with `flye`. Note that the `ont-hq` flag is used here, assuming you used Guppy5+ or Q20 (<5% error).
5. Visualize

## Installation
1. Create virtual environment and install all dependencies. Requires conda to be installed. See [here](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links) for instructions if needed:
```
# Create environment
conda create -n nanorgassm -y -c bioconda porechop filtlong minimap2 samtools flye bbmap git

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
2. Run NanOrgAssM to extract and assemble organelle reads matching your reference. Here we're going to extract chloroplast reads from garlic.
```
# Change path according to where you cloned NanOrgAssM on your machine.
# Don't forget to activate your nanorgassm virtual environment if not already done.
python nanorgassm.py \
  -r /path/to/my_reference.fasta \
  -i /path/to/my_reads.fastq.gz \
  -o /path/to/nanorgassm/output/folder \
  -t $(nproc)
  
```
`-r` is the reference sequence used to bait the reads
`-i` is the input reads. Fastq and fasta are supported, gzipped or not.
`-o` is the output folder.
`-t` is the number of threads used. `$(nproc)` means to use all availble threads (and should not be used on shared computing resources).

## Usage

