"""End-to-end pipeline tests with stub programs that mimic minimap2, filtlong, flye, ska and FastTree."""

import gzip
import json
import os
import stat
import sys
import textwrap
from pathlib import Path

import pytest

from bacon import BaconError
from bacon.pipeline import Settings, run
from bacon.seqio import read_records

STUBS = {
    # Reads whose name starts with "on" align; writes PAF lines. Samples named "bad*" crash.
    "minimap2": """
        import gzip, sys
        files = [a for a in sys.argv[1:] if a.endswith((".gz", ".fastq", ".fq", ".fasta"))][1:]
        if any("bad" in f for f in files):
            sys.exit("simulated crash")
        for f in files:
            with (gzip.open(f, "rt") if f.endswith(".gz") else open(f)) as fh:
                for i, line in enumerate(fh):
                    if i % 4 == 0 and line[1:].startswith("on"):
                        print(line[1:].split()[0] + "\\t100\\t0\\t100\\t+\\tref\\t1000\\t0\\t100\\t90\\t100\\t60")
    """,
    # Keeps reads at least --min_length long.
    "filtlong": """
        import gzip, sys
        args = sys.argv[1:]
        min_len = int(args[args.index("--min_length") + 1])
        with gzip.open(args[-1], "rt") as fh:
            lines = fh.read().splitlines()
        for i in range(0, len(lines), 4):
            if len(lines[i + 1]) >= min_len:
                print("\\n".join(lines[i:i + 4]))
    """,
    # Writes a one-contig circular assembly made of the first read.
    "flye": """
        import gzip, os, sys
        args = sys.argv[1:]
        out = args[args.index("--out-dir") + 1]
        os.mkdir(out)
        with gzip.open(args[args.index("--nano-hq") + 1], "rt") as fh:
            seq = fh.read().splitlines()[1]
        open(os.path.join(out, "assembly.fasta"), "w").write(">contig_1\\n" + seq + "\\n")
        open(os.path.join(out, "assembly_info.txt"), "w").write(
            "#seq_name\\tlength\\tcov.\\tcirc.\\trepeat\\tmult.\\talt_group\\tgraph_path\\n"
            f"contig_1\\t{len(seq)}\\t30\\tY\\tN\\t1\\t*\\t1\\n")
        open(os.path.join(out, "assembly_graph.gfa"), "w").write("H\\tVN:Z:1.0\\n")
    """,
    # ska build: nothing; ska align: an alignment of every genome in the build table (last base = sample index).
    "ska": """
        import os, sys
        args = sys.argv[1:]
        if args[0] == "build":
            prefix = args[args.index("-o") + 1]
            table = args[args.index("-f") + 1]
            open(prefix + ".skf", "w").write(open(table).read())
        else:
            out = args[args.index("-o") + 1]
            names = [line.split("\\t")[0] for line in open(args[-1]) if line.strip()]
            with open(out, "w") as fh:
                for i, name in enumerate(names):
                    seq = "" if os.environ.get("STUB_SKA_EMPTY") else f"ACGT{'ACGT'[i % 4]}"
                    fh.write(f">{name}\\n{seq}\\n")
    """,
    # sort: consume the SAM; index: nothing; consensus: the reference with its first 10 bases unknown (N).
    "samtools": """
        import sys
        from pathlib import Path
        args = sys.argv[1:]
        if args[0] == "sort":
            sys.stdin.read()
            open(args[args.index("-o") + 1], "w").write("BAM")
        elif args[0] == "consensus":
            out = Path(args[args.index("-o") + 1])
            ref = out.parents[2] / "reference.fasta"
            name, seq = ref.read_text().split("\\n", 1)
            seq = seq.replace("\\n", "")
            out.write_text(f"{name}\\n{'N' * 10}{seq[10:]}\\n")
    """,
    "Bandage": """
        import sys
        open(sys.argv[-1], "wb").write(b"PNG")
    """,
    "FastTree": """
        import sys
        names = [l[1:].split()[0] for l in open(sys.argv[-1]) if l.startswith(">")]
        print("(" + ",".join(f"{n}:0.{i + 1}" for i, n in enumerate(names)) + ");")
    """,
}


@pytest.fixture
def stubs(tmp_path, monkeypatch):
    """Put the stub programs first on PATH; every call is appended to calls.log."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    calls = tmp_path / "calls.log"
    for name, body in STUBS.items():
        path = bin_dir / name
        path.write_text(f"#!{sys.executable}\n"
                        f"import sys\n"
                        f"if sys.argv[1:] in (['--version'], ['-expert']):\n    sys.exit()\n"
                        f"open({str(calls)!r}, 'a').write({name!r} + '\\n')\n"
                        + textwrap.dedent(body))
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    return calls


def _reads(path: Path, reads: list[tuple[str, int]]) -> None:
    with gzip.open(path, "wt") as fh:
        for name, length in reads:
            fh.write(f"@{name}\n{'A' * length}\n+\n{'I' * length}\n")


@pytest.fixture
def dataset(tmp_path):
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\n" + "ACGT" * 250 + "\n")
    reads = tmp_path / "reads"
    reads.mkdir()
    for s in ("s1", "s2", "s3"):
        _reads(reads / f"{s}.fastq.gz", [("on1", 900), ("on2", 800), ("off1", 700)])
    _reads(reads / "none.fastq.gz", [("off1", 900)])  # Nothing matches the reference
    return ref, reads


def settings(ref: Path, reads: Path, out: Path, **kw) -> Settings:
    base = {"reference": ref, "input": reads, "output": out, "assembler": "flye", "snp_method": "ska", "threads": 2,
            "parallel": 2, "min_read_length": 100}
    return Settings(**{**base, **kw})


def calls(log: Path) -> list[str]:
    return log.read_text().split() if log.exists() else []


def summary(out: Path) -> dict[str, dict[str, str]]:
    lines = (out / "summary.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    return {row[0]: dict(zip(header, row)) for row in (line.split("\t") for line in lines[1:])}


def test_full_run(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    assert run(settings(ref, reads, out)) == 0
    rows = summary(out)
    assert rows["s1"]["Status"] == "ok"
    assert rows["s1"]["Raw_reads"] == "3" and rows["s1"]["Baited_reads"] == "2"
    assert rows["s1"]["Contigs"] == "1" and rows["s1"]["Circular_contigs"] == "1"
    assert rows["s1"]["Filtered_N50"] == "900"
    assert rows["none"]["Status"] == "failed (bait)"
    assert "no reads matched" in rows["none"]["Note"]
    assert [r.name for r in read_records(out / "3_assembled" / "all_assemblies" / "s1.fasta")] == ["s1_contig_1"]
    tree = (out / "4_compared" / "ska" / "tree.nwk").read_text()
    assert all(n in tree for n in ("s1", "s2", "s3", "Reference"))
    assert (out / "4_compared" / "ska" / "tree.svg").read_text().startswith("<svg")
    dist = (out / "4_compared" / "ska" / "snp_distances.tsv").read_text().splitlines()
    assert dist[0].split("\t") == ["snp-dists", "Reference", "s1", "s2", "s3"]
    info = json.loads((out / "run_info.json").read_text())
    assert info["samples"]["none"]["status"] == "failed (bait)"
    assert info["comparison"]["method"] == "ska"
    assert "flye" in info["tools"]
    assert "Baiting reads matching the reference" in (out / "bacon.log").read_text()  # INFO, not only warnings


def test_resume_skips_finished_steps(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    first = calls(stubs)
    stubs.unlink()
    run(settings(ref, reads, out))
    # Only the sample that failed at baiting is retried; nothing else runs again.
    assert calls(stubs) == ["minimap2"]
    assert len(first) > 5


def test_changed_parameter_reruns_that_step_and_later(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    stubs.unlink()
    run(settings(ref, reads, out, min_read_length=850))
    made = calls(stubs)
    assert made.count("filtlong") == 3 and made.count("flye") == 3 and "ska" in made
    assert made.count("minimap2") == 1  # Only the failed sample; baiting did not change
    assert summary(out)["s1"]["Filtered_reads"] == "1"


def test_changed_comparison_keeps_assemblies(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    stubs.unlink()
    run(settings(ref, reads, out, ska_min_freq=0.5))
    assert "flye" not in calls(stubs)
    assert (out / "4_compared" / "ska_0.5" / "tree.nwk").exists()
    assert (out / "4_compared" / "ska" / "tree.nwk").exists()  # The earlier comparison is kept


def test_redo(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    stubs.unlink()
    run(settings(ref, reads, out, redo="assemble"))
    made = calls(stubs)
    assert "filtlong" not in made and made.count("flye") == 3


def test_missing_output_triggers_rerun_of_that_sample(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    (out / "2_filtered" / "s2.fastq.gz").unlink()
    stubs.unlink()
    run(settings(ref, reads, out))
    made = calls(stubs)
    assert made.count("filtlong") == 1 and made.count("flye") == 1


def test_tool_crash_fails_only_that_sample(stubs, dataset, tmp_path):
    ref, reads = dataset
    _reads(reads / "bad.fastq.gz", [("on1", 900)])
    out = tmp_path / "out"
    assert run(settings(ref, reads, out)) == 0
    rows = summary(out)
    assert rows["bad"]["Status"] == "failed (bait)"
    assert "minimap2 failed" in rows["bad"]["Note"]
    assert "simulated crash" in (out / "logs" / "1_bait" / "bad.log").read_text()
    assert rows["s1"]["Status"] == "ok"


def test_fewer_than_three_assemblies_skips_comparison(stubs, dataset, tmp_path):
    ref, reads = dataset
    (reads / "s3.fastq.gz").unlink()
    out = tmp_path / "out"
    assert run(settings(ref, reads, out)) == 0
    assert not (out / "4_compared").exists()
    assert json.loads((out / "run_info.json").read_text())["comparison"] == {"skipped": "only 2 assemblies"}


def test_all_samples_failing_is_an_error(stubs, dataset, tmp_path):
    ref, reads = dataset
    for s in ("s1", "s2", "s3"):
        (reads / f"{s}.fastq.gz").unlink()
    with pytest.raises(BaconError, match="All samples failed at the bait step"):
        run(settings(ref, reads, tmp_path / "out"))


def test_missing_tool_is_reported(dataset, tmp_path, monkeypatch):
    ref, reads = dataset
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    with pytest.raises(BaconError, match="Required program.*minimap2.*conda install"):
        run(settings(ref, reads, tmp_path / "out"))


def test_low_depth_note(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out, snp_method="none"))
    rows = summary(out)
    assert rows["s1"]["Est_depth"] == "1.7"  # 1,700 bp of filtered reads over a 1 kb reference
    assert "low depth" in rows["s1"]["Note"]


def test_no_shared_snp_site_skips_the_tree(stubs, dataset, tmp_path, monkeypatch, caplog):
    ref, reads = dataset
    monkeypatch.setenv("STUB_SKA_EMPTY", "1")
    out = tmp_path / "out"
    assert run(settings(ref, reads, out)) == 0
    assert "try --ska-min-freq below 1" in caplog.text
    assert not (out / "4_compared" / "ska" / "tree.nwk").exists()
    assert (out / "4_compared" / "ska" / "snp_distances.tsv").exists()
    assert json.loads((out / "run_info.json").read_text())["comparison"]["tree"] is None


def test_failed_comparison_still_writes_the_summary(stubs, dataset, tmp_path):
    ref, reads = dataset
    (Path(os.environ["PATH"].split(os.pathsep)[0]) / "FastTree").write_text("#!/bin/sh\nexit 3\n")
    out = tmp_path / "out"
    with pytest.raises(BaconError, match="FastTree failed"):
        run(settings(ref, reads, out))
    assert summary(out)["s1"]["Status"] == "ok"
    assert "FastTree failed" in json.loads((out / "run_info.json").read_text())["comparison"]["failed"]


def test_templated_assembly(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    assert run(settings(ref, reads, out, assembler="samtools", snp_method="none")) == 0
    rows = summary(out)
    assert rows["s1"]["Assembly_length"] == "1000" and rows["s1"]["N_bases"] == "10"
    assert "10 N bases" in rows["s1"]["Note"]
    assert rows["s1"]["Circular_contigs"] == "NA"
    assert [r.name for r in read_records(out / "3_assembled" / "all_assemblies" / "s1.fasta")] == ["s1_ref"]
    # --template-gaps reference fills the leading N run from the reference.
    run(settings(ref, reads, out, assembler="samtools", snp_method="none", template_gaps="reference"))
    assert summary(out)["s1"]["N_bases"] == "0"


def test_adding_or_changing_a_sample_reruns_only_that_sample(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    stubs.unlink()
    _reads(reads / "s4.fastq.gz", [("on1", 900), ("on9", 850)])
    run(settings(ref, reads, out))
    made = calls(stubs)
    # s4 and the previously failed sample are baited; only s4 is filtered and assembled; all are compared.
    assert made.count("minimap2") == 2 and made.count("filtlong") == 1 and made.count("flye") == 1
    assert "ska" in made
    assert summary(out)["s4"]["Baited_reads"] == "2"
    stubs.unlink()
    os.utime(reads / "s2.fastq.gz", (1, 1))  # s2's reads changed (new modification time)
    run(settings(ref, reads, out))
    made = calls(stubs)
    assert made.count("filtlong") == 1 and made.count("flye") == 1


def test_add_genomes(stubs, dataset, tmp_path):
    ref, reads = dataset
    (reads / "s3.fastq.gz").unlink()  # Two samples + one added genome = three genomes
    public = tmp_path / "published.fa"
    public.write_text(">chr desc\nacgtacgt\n")
    out = tmp_path / "out"
    assert run(settings(ref, reads, out, add_genomes=[public])) == 0
    dist = (out / "4_compared" / "ska" / "snp_distances.tsv").read_text().splitlines()[0]
    assert dist.split("\t") == ["snp-dists", "Reference", "published", "s1", "s2"]
    assert (out / "4_compared" / "added_genomes" / "published.fasta").read_text().startswith(">published_chr\nACGT")
    clash = tmp_path / "s1.fasta"
    clash.write_text(">x\nACGT\n")
    with pytest.raises(BaconError, match="already used"):
        run(settings(ref, reads, out, add_genomes=[clash]))


def test_fasta_reads_are_filtered_without_filtlong(stubs, dataset, tmp_path):
    ref, reads = dataset
    for s in ("s1", "s2", "s3", "none"):
        (reads / f"{s}.fastq.gz").unlink()
    for s in ("f1", "f2", "f3"):
        (reads / f"{s}.fasta").write_text(">on1\n" + "A" * 900 + "\n>on2\n" + "A" * 50 + "\n")
    out = tmp_path / "out"
    assert run(settings(ref, reads, out)) == 0
    assert "filtlong" not in calls(stubs)
    assert summary(out)["f1"]["Filtered_reads"] == "1"
