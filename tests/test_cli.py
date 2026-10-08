import os
import re
from pathlib import Path

import pytest

import bacon
from bacon.cli import build_parser, main


def test_version_matches_pyproject():
    root = Path(__file__).parents[1]
    text = (root / "pyproject.toml").read_text()
    assert re.search(r'^version = "([^"]+)"', text, re.M).group(1) == bacon.__version__
    # The other files that give the version (see the release steps in docs/wiki/Development.md)
    v = bacon.__version__
    assert f"\nversion: {v}\n" in (root / "CITATION.cff").read_text()
    assert f"Nanopore reads (v{v})." in (root / "README.md").read_text()
    for doc in (root / "README.md", root / "docs" / "wiki" / "Installation.md"):  # Download and pip lines
        text = doc.read_text()
        assert set(re.findall(r"refs/tags/v([\d.]+)\.tar\.gz", text)) == {v}, doc
        assert set(re.findall(r"BACoN-([\d.]+)/example", text)) == {v}, doc
    for report in (root / "docs" / "reports").glob("*.html"):  # Published reports made with this version
        text = report.read_text()
        assert f"BACoN {v}" in text and not set(re.findall(r"BACoN (\d+\.\d+\.\d+)", text)) - {v}, report


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
    assert args.annotation is None


def test_annotation_option():
    args = build_parser().parse_args(["-r", "ref.gb", "-i", "in", "-o", "out", "--annotation", "ref.gff3"])
    assert args.annotation == Path("ref.gff3") and args.reference == Path("ref.gb")


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


def test_hdist_with_minimap2_is_ignored_with_a_warning(tmp_path, caplog, monkeypatch):
    import bacon.cli
    seen = []
    monkeypatch.setattr(bacon.cli, "run", lambda s: seen.append(s.hdist) or 0)
    base = ["-r", "r.fa", "-i", str(tmp_path), "-o", str(tmp_path / "o")]
    assert main([*base, "--hdist", "2"]) == 0
    assert "--hdist only applies to --baiting-method bbduk: ignored with minimap2" in caplog.text
    caplog.clear()
    assert main([*base, "-b", "bbduk", "--hdist", "2"]) == 0
    assert main([*base, "-b", "bbduk"]) == 0
    assert main(base) == 0
    assert "--hdist" not in caplog.text
    assert seen == [2, 2, 1, 1]


def test_an_interrupted_run_kills_the_programs_left_running(tmp_path, monkeypatch):
    import subprocess

    import bacon.cli
    from bacon import tools
    proc = subprocess.Popen(["sleep", "30"], start_new_session=True)
    with tools._RUNNING_LOCK:
        tools._RUNNING.add(proc)  # A program a second Ctrl-C kept the first one's cleanup from killing

    def interrupted(settings):
        raise KeyboardInterrupt

    monkeypatch.setattr(bacon.cli, "run", interrupted)
    try:
        assert main(["-r", "r.fa", "-i", str(tmp_path), "-o", str(tmp_path / "o")]) == 130
        assert proc.wait(timeout=5) == -9
    finally:
        proc.kill()
        with tools._RUNNING_LOCK:
            tools._RUNNING.discard(proc)


def test_after_the_run_ctrl_c_ends_the_program_without_a_traceback(tmp_path):
    import subprocess
    import sys
    # Run as a program (no arguments given to main): Python's exit, which waits for the threads, is not
    # interrupted by a KeyboardInterrupt
    code = ("import signal, sys, bacon.cli\n"
            "bacon.cli.run = lambda s: (_ for _ in ()).throw(KeyboardInterrupt)\n"
            f"sys.argv = ['bacon', '-r', 'r.fa', '-i', {str(tmp_path)!r}, '-o', {str(tmp_path / 'o')!r}]\n"
            "print(bacon.cli.main(), signal.getsignal(signal.SIGINT) is signal.SIG_DFL)\n")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout
    assert out.split() == ["130", "True"]
    # Called with arguments (as here), main leaves the handler of Ctrl-C as it was
    import signal
    before = signal.getsignal(signal.SIGINT)
    main(["-r", str(tmp_path / "missing.fa"), "-i", str(tmp_path), "-o", str(tmp_path / "o2"), "--snp-method", "none"])
    assert signal.getsignal(signal.SIGINT) is before


def test_bacon_error_exits_1(tmp_path, capsys):
    code = main(["-r", str(tmp_path / "missing.fa"), "-i", str(tmp_path), "-o", str(tmp_path / "o"),
                 "--snp-method", "none"])
    assert code == 1


def test_metadata_options():
    args = build_parser().parse_args(["-r", "r.fa", "-i", "in", "-o", "out", "--metadata", "m.tsv",
                                      "--color-by", "Group"])
    assert args.metadata == Path("m.tsv") and args.color_by == "Group"
    args = build_parser().parse_args(["-r", "r.fa", "-i", "in", "-o", "out"])
    assert args.metadata is None and args.color_by is None


def test_color_by_needs_metadata(capsys, tmp_path):
    from bacon.cli import main
    with pytest.raises(SystemExit):
        main(["-r", "r.fa", "-i", "in", "-o", "out", "--color-by", "group"])
    assert "--color-by needs --metadata" in capsys.readouterr().err
    # 'none' asks for nothing: no metadata needed (the run fails later, on the missing reference)
    assert main(["-r", str(tmp_path / "missing.fa"), "-i", str(tmp_path), "-o", str(tmp_path / "o"),
                 "--color-by", "None"]) == 1
    assert "--color-by needs" not in capsys.readouterr().err


def test_an_output_that_is_a_file_is_an_error_not_a_traceback(tmp_path, caplog):
    (tmp_path / "ref.fasta").write_text(">r\nACGT\n")
    (tmp_path / "out").write_text("x")
    assert main(["-r", str(tmp_path / "ref.fasta"), "-i", str(tmp_path / "ref.fasta"), "-o", str(tmp_path / "out"),
                 "--snp-method", "none"]) == 1
    assert f"The output folder {tmp_path / 'out'} is a file" in caplog.text


def test_changelog_has_this_version_and_the_citation_the_concept_doi():
    # A release renames "## Unreleased" to "## X.Y.Z (date)" (docs/wiki/Development.md); the release workflow
    # takes that section for the GitHub release
    root = Path(__file__).parents[1]
    assert re.search(rf"^## {re.escape(bacon.__version__)} \(\d{{4}}-\d{{2}}-\d{{2}}\)$",
                     (root / "CHANGELOG.md").read_text(), re.M)
    # "Cite this repository" at a tag gives the top-level DOI: the concept DOI, which resolves to the latest version
    citation = (root / "CITATION.cff").read_text()
    assert "\ndoi: 10.5281/zenodo.22970412\n" in citation and "Concept DOI" in citation


def test_color_by_help_gives_the_limit_of_values(capsys):
    with pytest.raises(SystemExit):
        build_parser().parse_args(["--help"])
    text = " ".join(capsys.readouterr().out.split())
    assert "at most 48 distinct values (each gets a colour and a shape)" in text
    # Not only the tree and the heatmap: every figure and table that shows the genomes
    assert text.count("colours every genome in the report (tree, heatmap, bar charts, tables)") == 2


@pytest.mark.skipif(os.geteuid() == 0, reason="root reads anything")
def test_a_path_that_cannot_be_read_is_an_error_not_a_traceback(tmp_path, caplog, monkeypatch):
    import bacon.pipeline
    monkeypatch.setattr(bacon.pipeline, "require", lambda names: {})  # The programs are not needed here
    (tmp_path / "ref.fasta").write_text(">r\nACGT\n")
    locked = tmp_path / "locked"
    locked.mkdir()
    locked.chmod(0)
    ref = str(tmp_path / "ref.fasta")
    try:
        for argv, path in ((["-i", str(locked / "reads"), "-o", str(tmp_path / "out")], locked / "reads"),
                           (["-i", ref, "-o", str(locked / "out")], locked / "out")):
            caplog.clear()
            assert main(["-r", ref, *argv, "--snp-method", "none"]) == 1
            assert [r.getMessage() for r in caplog.records] == [f"{path}: Permission denied"]
        # In a started run, the error is in bacon.log too
        caplog.clear()
        out = tmp_path / "out"
        assert main(["-r", str(locked / "ref.fasta"), "-i", ref, "-o", str(out), "--snp-method", "none"]) == 1
        assert f"{locked / 'ref.fasta'}: Permission denied" in caplog.text
        assert f"[ERROR] {locked / 'ref.fasta'}: Permission denied" in (out / "bacon.log").read_text()
    finally:
        locked.chmod(0o700)
