# Installation

BACoN is a Python package (standard library only) that runs command-line programs, which are installed with
conda.

## bioconda (recommended)

```bash
conda create -n bacon -c conda-forge -c bioconda bacon-nanopore
conda activate bacon
bacon --version
```

The package is named `bacon-nanopore` because conda-forge already has an unrelated package called `bacon`; the
command is `bacon`. It installs BACoN with every program below except Bandage.

## From the source code

For development, or the latest code:

```bash
git clone https://github.com/duceppemo/BACoN
cd BACoN
conda env create -f environment.yml
conda activate BACoN
pip install .
bacon --version
```

[`environment.yml`](https://github.com/duceppemo/BACoN/blob/main/environment.yml) lists the programs and their
minimum versions:

| Program | Used for | Minimum |
|---|---|---|
| minimap2 | baiting; templated assembly | 2.24 |
| samtools | templated assembly (`samtools consensus -X r10.4_sup`) | 1.21 |
| Filtlong | read filtering | 0.2.1 |
| Flye | de novo assembly (`-a flye`) | 2.9.5 |
| myloasm | de novo assembly (`-a myloasm`) | 0.7 |
| SKA2 | SNPs (`--snp-method ska`, default) | 0.5 |
| Parsnp, HarvestTools | SNPs (`--snp-method parsnp`) | 2.1.2 (with 2.1.1, every SNP is one base off) |
| FastTree | tree (default) | 2.1.11 |
| IQ-TREE | tree (`--tree iqtree`) | 2.2 |
| BBMap (BBDuk) | baiting (`-b bbduk`) | 39 |
| Bandage | optional: pictures of the assembly graphs | |
| python-isal | optional: faster reading of gzipped reads | |

BACoN checks at start-up that the programs needed by the chosen options are on the `PATH`, and lists the
missing ones with the command to install them.

## Checking the installation

The example is in the repository (`example/`), not in the conda package. From a clone, or after downloading the
example folder of the release:

```bash
curl -sL https://github.com/duceppemo/BACoN/archive/refs/tags/v0.3.8.tar.gz | tar -xz --strip-components=1 BACoN-0.3.8/example
bash example/run_example.sh
```

generates a small simulated dataset (an annotated 30 kb circular reference, four samples with known SNPs and
their metadata), runs BACoN with the default settings and checks the results against the truth: every pairwise
SNP distance, the region, gene, context and effect of every SNP, and the report's annotation and metadata
([Example](Example)); it prints one `OK` line per check. It takes a few seconds.

## pip only

If the programs are already installed (for example in an HPC module system), BACoN itself installs with
`pip install https://github.com/duceppemo/BACoN/archive/refs/tags/v0.3.8.tar.gz` (or
`pip install git+https://github.com/duceppemo/BACoN` for the latest code).

## Troubleshooting

- **The environment takes long to solve, or mamba reports "nothing provides ..." for packages that exist**:
  some mamba 2.x versions misreport conflicts; `conda env create` (with the libmamba solver) solves it.
- **`samtools consensus: unrecognised option -X`**: samtools is older than 1.17 (which added `-X`); update
  it. BACoN is tested with samtools 1.21 and later.
- **Upgrading from BACoN 0.2**: create a new environment; the old `requirements.txt` environment pinned
  programs that are no longer used (Porechop, Shasta, Rebaler, Snippy, PhaME, RAxML, ete3).
