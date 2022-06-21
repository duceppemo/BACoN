# NanOrgAssM
NanOrgAssM stands for "Nanopore Organelle Assembly Method". As its name suggests, it's aimed at assembling organelle genomes from Nanopore genome skimming. Hopefully, targeted sequencing with a suitable organelle reference sequencing has been used to enrich organelle sequences in the reads (I like to call that "powerskimming").

## TODO

## Information
NanoOrgAssM follows these steps:
1. Basecall and demultiplex the fast5 with `guppy` in SUP mode to generate the fastq files.
2. Remove adapters with `porechop` and filter low quality reads with `filtlong`.
3. Extract organelle reads with `minimap2` using a provided reference sequence (ideally the same used for the targeted sequencing).
4. Assemble the extracted reads with `flye`.
5. Visualize

## Usage

## Example
