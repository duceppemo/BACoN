import re
from pathlib import Path

import pytest

import bacon
from bacon.cli import build_parser, main


def test_version_matches_pyproject():
    text = (Path(__file__).parents[1] / "pyproject.toml").read_text()
    assert re.search(r'^version = "([^"]+)"', text, re.M).group(1) == bacon.__version__


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert capsys.readouterr().out.strip() == f"BACoN {bacon.__version__}"


def test_defaults():
    args = build_parser().parse_args(["-r", "ref.fa", "-i", "in", "-o", "out"])
    assert (args.baiting_method, args.assembly_method, args.snp_method, args.tree) == (
        "minimap2", "samtools", "ska", "fasttree")
    assert args.kmer_size == 31 and args.min_read_length == 500 and args.ska_min_freq == 1.0


def test_legacy_snp_flag():
    args = build_parser().parse_args(["-r", "r", "-i", "i", "-o", "o", "-snp", "parsnp"])
    assert args.snp_method == "parsnp"


@pytest.mark.parametrize("argv, message", [
    (["-r", "r", "-o", "o"], "exactly one of -i/--input or --sample-sheet"),
    (["-r", "r", "-i", "i", "--sample-sheet", "s", "-o", "o"], "exactly one of"),
    (["-r", "r", "-i", "i", "-o", "o", "-k", "99"], "at most 31"),
    (["-r", "r", "-i", "i", "-o", "o", "--keep-percent", "0"], "above 0"),
    (["-r", "r", "-i", "i", "-o", "o", "--ska-min-freq", "1.5"], "at most 1"),
    (["-r", "r", "-i", "i", "-o", "o", "-t", "0"], "at least 1"),
    (["-r", "r", "-i", "i", "-o", "o", "--keep-bam", "-b", "bbduk"], "--keep-bam only applies"),
    (["-r", "r", "-i", "i", "-o", "o", "-a", "shasta"], "invalid choice"),
    (["-r", "r", "-i", "i", "-o", "o", "--snp-method", "snippy"], "invalid choice"),
])
def test_argument_errors(argv, message, capsys):
    with pytest.raises(SystemExit) as exc:
        main(argv)
    assert exc.value.code == 2
    assert message in capsys.readouterr().err


def test_bacon_error_exits_1(tmp_path, capsys):
    code = main(["-r", str(tmp_path / "missing.fa"), "-i", str(tmp_path), "-o", str(tmp_path / "o"),
                 "--snp-method", "none"])
    assert code == 1
