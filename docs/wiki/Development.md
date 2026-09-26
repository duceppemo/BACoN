# Development

## Setup

```bash
git clone https://github.com/duceppemo/BACoN && cd BACoN
conda env create -f environment.yml && conda activate BACoN
pip install -e ".[test]"
pre-commit install
```

## Code

| Module | Contents |
|---|---|
| `bacon/cli.py` | command line |
| `bacon/pipeline.py` | the steps, checkpoints, summary and provenance |
| `bacon/samples.py` | finding samples (file, folder, sample sheet) |
| `bacon/steps.py` | baiting, filtering, assembly |
| `bacon/compare.py` | SKA2, Parsnp, distances, trees |
| `bacon/newick.py` | Newick parsing, midpoint rooting, SVG drawing |
| `bacon/report.py` | the HTML report (`python -m bacon.report OUTPUT` rebuilds it) |
| `bacon/multiqc.py` | the MultiQC custom-content files |
| `bacon/seqio.py` | fasta/fastq reading and writing |
| `bacon/tools.py` | running external programs and logging their output |

BACoN uses the Python standard library only (python-isal is used when installed). External programs are run
with `bacon.tools.run`, which logs the command and the program's messages and raises an error with the end of
the log when a program fails.

## Tests

```bash
ruff check .
pytest --cov=bacon
```

The tests replace the external programs by small stub scripts (`tests/test_pipeline.py`), so they run in
seconds, without conda, and test the pipeline logic: resuming, parameter changes, failures, provenance. The
real programs are exercised by the example and the validation:

```bash
bash example/run_example.sh                   # seconds; checks SNP distances against the truth
bash validation/run_validation.sh /tmp/val    # a few minutes; scores every assembler and SNP method
```

Continuous integration (`.github/workflows/ci.yml`) runs ruff, the tests on Python 3.10 to 3.13 (and macOS),
and the example in the conda environment.

## Adding or replacing a program

A program is added or replaces another only if it is maintained, installs with the others from bioconda,
and does at least as well in the validation. Add it to `environment.yml`, `recipe/meta.yaml` and
`bacon/tools.py` (`PACKAGES`), with a stub in `tests/test_pipeline.py`, then run the validation and record the
results in a new dated folder under `validation/results/` (never rewrite an old record).

## Documentation

The wiki is maintained in `docs/wiki/` and published to the GitHub wiki by `.github/workflows/wiki.yml` on
every push to `main` that changes it. Do not edit the wiki on GitHub: the next publication overwrites it.

## Releases

1. Update `version` in `pyproject.toml`, `__version__` in `bacon/__init__.py`, `CITATION.cff` (version and
   date) and `CHANGELOG.md`.
2. Commit, tag `vX.Y.Z` and push the tag: `.github/workflows/release.yml` checks the versions, builds the
   package and creates the GitHub release with the changelog section.
3. Zenodo archives the release and mints a version DOI: add it to `CITATION.cff` (`doi` and `identifiers`).
   The README badge uses the concept DOI (10.5281/zenodo.22970412), which always points to the latest
   version.
4. Update `recipe/meta.yaml` (version, sha256 of the tag's tarball) and the bioconda-recipes pull request.
