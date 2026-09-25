# Contributing

Thanks for helping improve BACoN. Bug reports, validation results and documentation fixes are as welcome as code.

## Reporting a bug

Open an [issue](https://github.com/duceppemo/BACoN/issues/new/choose) with the exact command, `bacon.log` and
`run_info.json` from the output folder (they record the versions of BACoN and of every program it ran), and the
log of the step that failed (`logs/<step>/<sample>.log`).

## Changing code

```bash
git clone https://github.com/duceppemo/BACoN && cd BACoN
conda env create -f environment.yml && conda activate BACoN
pip install -e ".[test]"
pre-commit install
```

- Keep `ruff check .` and `pytest` green; add a test for every fix or feature. The tests replace the external
  programs with stubs (`tests/test_pipeline.py`), so they run in seconds without conda.
- BACoN itself uses the Python standard library only; external programs are called on the command line.
- A new or replacement program must be maintained, installable from bioconda together with the others, and
  validated: run `validation/run_validation.sh` and compare its scores with the current tools.
- User-visible changes go in `CHANGELOG.md`; options and outputs are documented in `docs/wiki/`.
- Open the pull request against `main`.
