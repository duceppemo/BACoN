<p align="center">
  <img src="docs/images/bacon_logo.png" alt="BACoN: Bait, Assemble and Compare Nanopore" width="600">
</p>

<p align="center">
  <a href="https://github.com/duceppemo/BACoN/actions/workflows/ci.yml"><img src="https://github.com/duceppemo/BACoN/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://codecov.io/gh/duceppemo/BACoN"><img src="https://codecov.io/gh/duceppemo/BACoN/graph/badge.svg" alt="Coverage"></a>
  <a href="https://github.com/duceppemo/BACoN/releases/latest"><img src="https://img.shields.io/github/v/release/duceppemo/BACoN?label=release&cacheSeconds=3600" alt="Latest release"></a>
  <img src="https://img.shields.io/badge/python-3.10%2B-blue" alt="Python 3.10+">
  <a href="LICENSE"><img src="https://img.shields.io/github/license/duceppemo/BACoN" alt="License: MIT"></a>
  <a href="https://github.com/duceppemo/BACoN/wiki"><img src="https://img.shields.io/badge/docs-wiki-informational" alt="Documentation"></a>
  <a href="https://doi.org/10.5281/zenodo.22970412"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.22970412.svg" alt="DOI"></a>
</p>

BACoN extracts the Nanopore reads that match a reference sequence (typically an organelle genome) from many
samples, assembles them, and compares the samples: SNP distances and a tree. It was designed for genome
skimming, and works best with targeted sequencing that enriches the region of interest.

```
reads ──► 1. bait (minimap2) ──► 2. filter (Filtlong) ──► 3. assemble ──────────────► 4. compare
                                                          templated: samtools consensus    SNPs: SKA2 or Parsnp
                                                          de novo:   Flye or myloasm       tree: FastTree or IQ-TREE
```

## Quick start

```bash
git clone https://github.com/duceppemo/BACoN && cd BACoN
conda env create -f environment.yml && conda activate BACoN
pip install .

bacon -r chloroplast.fasta -i reads/ -o results/
```

`reads/` holds one fastq file per sample, or one folder per sample (such as MinKNOW's `fastq_pass/barcode01/`).
Open `results/report.html` for an overview: samples with their read counts, depth and assembly, the tree, a
heatmap of the SNP distances, and a methods paragraph. The same results are in `summary.tsv`,
`4_compared/ska/snp_distances.tsv` and `4_compared/ska/tree.nwk`, and MultiQC picks up BACoN's
`*_mqc.json` files. Rerunning the same command resumes where it stopped; changing a parameter reruns only the
steps it affects.

To check an installation, run the bundled example (simulated reads with known SNPs; a few seconds):

```bash
bash example/run_example.sh
```

## Documentation

Everything else is in the [wiki](https://github.com/duceppemo/BACoN/wiki), whose sources are maintained in
[`docs/wiki`](docs/wiki):

- [Installation](https://github.com/duceppemo/BACoN/wiki/Installation) — conda, pip, troubleshooting
- [Usage](https://github.com/duceppemo/BACoN/wiki/Usage) — inputs, options, resuming, performance
- [Methods](https://github.com/duceppemo/BACoN/wiki/Methods) — what each step does, and which assembler or SNP method to choose
- [Tutorial](https://github.com/duceppemo/BACoN/wiki/Tutorial) — 28 potato cultivars from a public project
- [Outputs](https://github.com/duceppemo/BACoN/wiki/Outputs) — files and columns
- [Validation](https://github.com/duceppemo/BACoN/wiki/Validation) — how the methods were chosen and tested
- [FAQ](https://github.com/duceppemo/BACoN/wiki/FAQ) — troubleshooting
- [Development](https://github.com/duceppemo/BACoN/wiki/Development) — tests, CI, releases

## Citation

If BACoN is useful in your work, please cite it (see [`CITATION.cff`](CITATION.cff)):

> Duceppe, M.-O. (2026). BACoN: Bait, Assemble and Compare Nanopore reads (v0.3.0). Zenodo. https://doi.org/10.5281/zenodo.22970412

and the programs it runs: minimap2, Filtlong, samtools or Flye or myloasm, SKA2 or Parsnp, and FastTree or
IQ-TREE ([references](https://github.com/duceppemo/BACoN/wiki/Methods#references)).

## License

[MIT](LICENSE)
