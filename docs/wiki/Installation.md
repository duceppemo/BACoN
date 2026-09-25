# Installation

BACoN is a Python package (standard library only) that runs command-line programs, which are installed with
conda.

## conda (recommended)

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
| Parsnp, HarvestTools | SNPs (`--snp-method parsnp`) | 2.0 |
| FastTree | tree (default) | 2.1.11 |
| IQ-TREE | tree (`--tree iqtree`) | 2.2 |
| BBMap (BBDuk) | baiting (`-b bbduk`) | 39 |
| Bandage | optional: pictures of the assembly graphs | |
| python-isal | optional: faster reading of gzipped reads | |

BACoN checks at start-up that the programs needed by the chosen options are on the `PATH`, and lists the
missing ones with the command to install them.

## Checking the installation

```bash
bash example/run_example.sh
```

generates a small simulated dataset (a 30 kb circular reference and four samples with known SNPs), runs BACoN
with the default settings and checks every pairwise SNP distance against the truth. It takes a few
seconds.

## pip only

If the programs are already installed (for example in an HPC module system), BACoN itself installs with
`pip install https://github.com/duceppemo/BACoN/archive/refs/tags/v0.3.0.tar.gz`.

## Troubleshooting

- **The environment takes long to solve, or mamba reports "nothing provides ..." for packages that exist**:
  some mamba 2.x versions misreport conflicts; `conda env create` (with the libmamba solver) solves it.
- **`samtools consensus: unrecognised option -X`**: samtools is older than 1.17 (which added `-X`); update
  it. BACoN is tested with samtools 1.21 and later.
- **Upgrading from BACoN 0.2**: create a new environment; the old `requirements.txt` environment pinned
  programs that are no longer used (Porechop, Shasta, Rebaler, Snippy, PhaME, RAxML, ete3).
