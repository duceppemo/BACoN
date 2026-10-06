"""End-to-end pipeline tests with stub programs that mimic minimap2, filtlong, flye, ska and FastTree."""

import gzip
import json
import os
import stat
import sys
import textwrap
from pathlib import Path

import pytest

from bacon import BaconError, __version__
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
        elif args[0] == "map" and os.environ.get("STUB_SKA_MAP_FAIL"):
            open(args[args.index("-o") + 1], "w").write("partial")
            sys.exit("map failed")
        elif args[0] == "map":  # VCF: one variant, carried by every genome but the reference
            names = [line.split("\\t")[0] for line in open(args[2]) if line.strip()]
            header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *names]
            record = ["ref", "5", ".", "A", "C", ".", ".", ".", "GT"]
            record += ["0" if n == "Reference" else "1" for n in names]
            with open(args[args.index("-o") + 1], "w") as fh:
                fh.write("##fileformat=VCFv4.4\\n##contig=<ID=ref>\\n")
                fh.write("\\t".join(header) + "\\n" + "\\t".join(record) + "\\n")
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
        if args[0] == "sort":  # The "BAM" holds the aligner's output as is
            open(args[args.index("-o") + 1], "w").write(sys.stdin.read())
        elif args[0] == "view":
            sys.stdout.write(open(args[-1]).read())
        elif args[0] == "consensus":
            out = Path(args[args.index("-o") + 1])
            ref = out.parents[2] / "reference.fasta"
            name, seq = ref.read_text().split("\\n", 1)
            seq = seq.replace("\\n", "")
            out.write_text(f"{name}\\n{'N' * 10}{seq[10:]}\\n")
    """,
    # Keeps reads whose name starts with "on"; reports BBDuk's counts; STUB_BBDUK_OOM: runs out of memory.
    "bbduk.sh": """
        import gzip, os, sys
        opts = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
        if os.environ.get("STUB_BBDUK_OOM"):
            sys.stderr.write("Allocating kmer table: Terminating due to java.lang.OutOfMemoryError: Java heap space\\n")
            sys.exit(3)
        with gzip.open(opts["in"], "rt") as fh:
            lines = fh.read().splitlines()
        kept = [lines[i:i + 4] for i in range(0, len(lines), 4) if lines[i][1:].startswith("on")]
        with gzip.open(opts["outm"], "wt") as out:
            out.write("".join("\\n".join(r) + "\\n" for r in kept))
        total = sum(len(lines[i + 1]) for i in range(0, len(lines), 4))
        sys.stderr.write(f"Input:                  \\t{len(lines) // 4} reads \\t\\t{total} bases.\\n")
        sys.stderr.write(f"Contaminants:           \\t{len(kept)} reads (x%) \\t{sum(len(r[1]) for r in kept)} bases (x%)\\n")
        sys.stderr.write("hdist=" + opts["hdist"] + "\\n")
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
    vcf = (out / "4_compared" / "ska" / "snps.vcf").read_text().splitlines()
    assert vcf[0] == "##fileformat=VCFv4.2" and "##contig=<ID=ref,length=1000>" in vcf
    assert vcf[-2].split("\t")[9:] == ["s1", "s2", "s3"] and vcf[-1].split("\t")[9:] == ["1", "1", "1"]
    assert info["comparison"]["vcf"].endswith("snps.vcf")
    import hashlib
    assert info["reference"]["md5"] == hashlib.md5(ref.read_bytes()).hexdigest()  # The file given, not the copy
    assert info["reference"]["md5"] != hashlib.md5((out / "reference.fasta").read_bytes()).hexdigest()
    report = (out / "report.html").read_text()
    assert all(name in report for name in ("s1", "s2", "s3", "none", "SKA2", "SNP distances", "4 genomes, 4 distinct"))
    for name in ("bacon_samples_mqc.json", "bacon_reads_mqc.json", "bacon_distances_mqc.json"):
        json.loads((out / name).read_text())
    assert set(json.loads((out / "bacon_samples_mqc.json").read_text())["data"]) == {"s1", "s2", "s3", "none"}


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


def test_report_failure_does_not_fail_the_run(stubs, dataset, tmp_path, monkeypatch, caplog):
    import bacon.pipeline

    def broken(output):
        raise RuntimeError("boom")

    monkeypatch.setattr(bacon.pipeline, "write_report", broken)
    ref, reads = dataset
    out = tmp_path / "out"
    assert run(settings(ref, reads, out, snp_method="none")) == 0
    assert "Could not write the HTML report: boom" in caplog.text
    assert (out / "summary.tsv").exists() and not (out / "bacon_distances_mqc.json").exists()


def test_tree_keeps_names_the_tree_programs_would_change(stubs, dataset, tmp_path):
    ref, reads = dataset
    (reads / "s1.fastq.gz").rename(reads / "a+b.fastq.gz")
    out = tmp_path / "out"
    assert run(settings(ref, reads, out)) == 0
    tree = (out / "4_compared" / "ska" / "tree.nwk").read_text()
    assert "a+b" in tree and "g0000" not in tree
    assert "a+b" in (out / "4_compared" / "ska" / "tree.svg").read_text()
    assert not (out / "4_compared" / "ska" / "tree_input.fasta").exists()


def test_notes_with_tabs_or_line_breaks_keep_the_summary_rectangular(tmp_path):
    from bacon.pipeline import SampleState, summary_row
    from bacon.samples import Sample
    st = SampleState(Sample("x", [tmp_path / "x.fq"]), failed="failed (bait)",
                     notes=["error:\tbad\ninput", "second"])
    assert summary_row(st)["Note"] == "error: bad input; second"


def test_vcf_failure_is_a_warning_and_a_rerun_makes_the_vcf(stubs, dataset, tmp_path, monkeypatch, caplog):
    ref, reads = dataset
    out = tmp_path / "out"
    ska_dir = out / "4_compared" / "ska"
    monkeypatch.setenv("STUB_SKA_MAP_FAIL", "1")
    assert run(settings(ref, reads, out)) == 0
    assert "Could not write the VCF" in caplog.text
    assert not (ska_dir / "snps.vcf").exists() and not (ska_dir / "snps.raw.vcf").exists()
    assert json.loads((out / "run_info.json").read_text())["comparison"]["vcf"] is None
    monkeypatch.delenv("STUB_SKA_MAP_FAIL")
    stubs.unlink()
    assert run(settings(ref, reads, out)) == 0  # The comparison is reused; only the VCF is made
    assert calls(stubs).count("ska") == 1 and "FastTree" not in calls(stubs)
    assert (ska_dir / "snps.vcf").exists()
    assert json.loads((out / "run_info.json").read_text())["comparison"]["vcf"].endswith("snps.vcf")


def test_comparison_from_before_the_vcf_gets_one(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    ska_dir = out / "4_compared" / "ska"
    (ska_dir / "snps.vcf").unlink()
    checkpoint = out / ".checkpoints" / "compare.json"
    data = json.loads(checkpoint.read_text())
    del data["results"]["vcf"]  # As written by BACoN 0.3.1, which kept no ska_reference.fasta either
    checkpoint.write_text(json.dumps(data))
    (ska_dir / "ska_reference.fasta").unlink()
    run(settings(ref, reads, out))
    assert (ska_dir / "snps.vcf").exists() and (ska_dir / "ska_reference.fasta").exists()
    assert "vcf" in json.loads(checkpoint.read_text())["results"]


def test_vcf_from_before_0_3_3_is_rewritten_on_resume(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    vcf = out / "4_compared" / "ska" / "snps.vcf"
    current = vcf.read_text()
    vcf.write_text(current.replace(f"##source=BACoN {__version__}", "##source=BACoN 0.3.2"))
    stubs.unlink()
    run(settings(ref, reads, out))
    assert vcf.read_text() == current  # Rewritten from the comparison's files, without redoing the comparison
    assert calls(stubs).count("ska") == 1 and "FastTree" not in calls(stubs)
    stubs.unlink()
    run(settings(ref, reads, out))
    assert "ska" not in calls(stubs)  # A current VCF is kept


def test_reference_and_added_genomes_are_upper_case_acgtn(stubs, dataset, tmp_path):
    ref, reads = dataset
    soft = tmp_path / "soft.fasta"
    soft.write_text(">ref masked\n" + "acgt" * 10 + "RYKM" + "ACGT" * 240 + "\n")
    public = tmp_path / "pub.fasta"
    public.write_text(">p\nacgtnRYacgt\n")
    out = tmp_path / "out"
    assert run(settings(soft, reads, out, add_genomes=[public])) == 0
    copy = "".join((out / "reference.fasta").read_text().split("\n")[1:])
    assert copy == "ACGT" * 10 + "NNNN" + "ACGT" * 240
    assert (out / "4_compared" / "added_genomes" / "pub.fasta").read_text() == ">pub_p\nACGTNNNACGT\n"


def test_added_genomes_do_not_redo_the_comparison(stubs, dataset, tmp_path):
    ref, reads = dataset
    public = tmp_path / "pub.fasta"
    public.write_text(">p\nACGTACGT\n")
    out = tmp_path / "out"
    run(settings(ref, reads, out, add_genomes=[public]))
    stubs.unlink()
    import time
    time.sleep(1.1)  # A new modification time for anything rewritten
    run(settings(ref, reads, out, add_genomes=[public]))
    assert "ska" not in calls(stubs) and "FastTree" not in calls(stubs)
    public.write_text(">p\nACGTACGA\n")  # A changed genome: compared again
    run(settings(ref, reads, out, add_genomes=[public]))
    assert "ska" in calls(stubs)


def test_unexpected_error_fails_one_sample_not_the_run(stubs, dataset, tmp_path, monkeypatch):
    import bacon.steps
    real = bacon.steps.filter_filtlong

    def flaky(name, *args, **kwargs):
        if name == "s2":
            raise UnicodeEncodeError("ascii", "x", 0, 1, "boom")
        return real(name, *args, **kwargs)

    monkeypatch.setattr(bacon.steps, "filter_filtlong", flaky)
    out = tmp_path / "out"
    assert run(settings(*dataset, out, snp_method="none")) == 0
    rows = summary(out)
    assert rows["s2"]["Status"] == "failed (filter)" and "unexpected error: UnicodeEncodeError" in rows["s2"]["Note"]
    sample_log = (out / "logs" / "2_filter" / "s2.log").read_text()
    assert "Traceback" in sample_log and "flaky" in sample_log  # Where the error happened, for a bug report
    assert rows["s1"]["Status"] == "ok"


def test_unexpected_error_in_the_comparison_keeps_the_summary(stubs, dataset, tmp_path, monkeypatch):
    import bacon.compare

    def broken(*args, **kwargs):
        raise FileNotFoundError("missing.contree")

    monkeypatch.setattr(bacon.compare, "run_ska", broken)
    out = tmp_path / "out"
    with pytest.raises(BaconError, match="comparison failed unexpectedly"):
        run(settings(*dataset, out))
    assert summary(out)["s1"]["Status"] == "ok"
    assert "FileNotFoundError" in json.loads((out / "run_info.json").read_text())["comparison"]["failed"]


def test_templated_assembly_removes_graphs_of_an_earlier_de_novo_assembly(stubs, dataset, tmp_path):
    out = tmp_path / "out"
    run(settings(*dataset, out, snp_method="none"))
    graphs = out / "3_assembled" / "assembly_graphs"
    assert (graphs / "s1.gfa").exists()
    run(settings(*dataset, out, snp_method="none", assembler="samtools"))
    assert not list(graphs.glob("s1.*"))


def test_an_interruption_loses_only_the_samples_still_running(stubs, dataset, tmp_path, monkeypatch):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out, snp_method="none"))
    _reads(reads / "s4.fastq.gz", [("on1", 900), ("on2", 850)])
    import bacon.steps
    real = bacon.steps.assemble_flye

    def interrupted(name, *args, **kwargs):
        if name == "s4":
            raise KeyboardInterrupt
        return real(name, *args, **kwargs)

    monkeypatch.setattr(bacon.steps, "assemble_flye", interrupted)
    with pytest.raises(KeyboardInterrupt):
        run(settings(ref, reads, out, snp_method="none"))
    monkeypatch.setattr(bacon.steps, "assemble_flye", real)
    stubs.unlink()
    run(settings(ref, reads, out, snp_method="none"))
    assert calls(stubs).count("flye") == 1  # Only s4: the assemblies of s1-s3 were kept


def test_bbduk_counts_hdist_and_out_of_memory(stubs, dataset, tmp_path, monkeypatch):
    ref, reads = dataset
    out = tmp_path / "out"
    assert run(settings(ref, reads, out, baiting="bbduk", hdist=1, snp_method="none")) == 0
    rows = summary(out)
    assert rows["s1"]["Raw_reads"] == "3" and rows["s1"]["Baited_reads"] == "2"
    assert "hdist=1" in (out / "logs" / "1_bait" / "s1.log").read_text()
    monkeypatch.setenv("STUB_BBDUK_OOM", "1")
    with pytest.raises(BaconError, match="All samples failed"):
        run(settings(ref, reads, tmp_path / "oom", baiting="bbduk", hdist=2, memory_gb=4, snp_method="none"))
    note = summary(tmp_path / "oom")["s1"]["Note"]
    assert "ran out of memory with 2 GB" in note and "--hdist 1" in note


def test_added_genome_cannot_take_the_name_of_a_failed_sample(stubs, dataset, tmp_path):
    ref, reads = dataset
    clash = tmp_path / "none.fasta"  # Sample "none" fails at baiting
    clash.write_text(">x\nACGT\n")
    with pytest.raises(BaconError, match="already used"):
        run(settings(ref, reads, tmp_path / "out", add_genomes=[clash]))


def test_interruption_kills_the_programs_still_running(tmp_path, caplog):
    import time

    from bacon.pipeline import SampleState, _run_parallel
    from bacon.samples import Sample
    from bacon.steps import StepResult
    from bacon.tools import run as run_tool

    pidfile = tmp_path / "child.pid"

    def fn(st, threads):
        if st.sample.name == "slow":  # A program that starts a program of its own, like Flye
            run_tool(["sh", "-c", f"sleep 30 & echo $! > {pidfile}.tmp; mv {pidfile}.tmp {pidfile}; wait"],
                     tmp_path / "slow.log")
        else:  # "fast" ends (and the interruption comes) once the slow program and its child are running
            for _ in range(100):
                if pidfile.exists():
                    break
                time.sleep(0.05)
        return StepResult(None)

    def on_done(name, res):
        if name == "fast":
            raise KeyboardInterrupt

    states = [SampleState(Sample(n, [tmp_path / f"{n}.fq"])) for n in ("slow", "fast")]
    start = time.monotonic()
    with pytest.raises(KeyboardInterrupt):
        _run_parallel(states, fn, settings(tmp_path, tmp_path, tmp_path, parallel=2), "bait", tmp_path,
                      on_done=on_done)
    assert time.monotonic() - start < 10  # The program was killed
    assert "slow:" not in caplog.text  # Killed by the interruption: not reported as a failure of the sample
    assert pidfile.exists()
    if pidfile.exists():  # ... and so was the program it started
        import os
        pid = int(pidfile.read_text())
        time.sleep(0.5)
        try:
            os.kill(pid, 0)
            alive = Path(f"/proc/{pid}/stat").read_text().split()[2] != "Z"
        except (ProcessLookupError, FileNotFoundError):
            alive = False
        assert not alive


def test_two_runs_cannot_share_an_output_folder(stubs, dataset, tmp_path):
    import fcntl
    out = tmp_path / "out"
    out.mkdir()
    with open(out / ".bacon.lock", "w") as held:
        fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)  # Another run holds it
        with pytest.raises(BaconError, match="Another BACoN run is using"):
            run(settings(*dataset, out))
    assert run(settings(*dataset, out, snp_method="none")) == 0  # Released: this run goes ahead


def test_keep_bam(stubs, dataset, tmp_path):
    out = tmp_path / "out"
    assert run(settings(*dataset, out, keep_bam=True, snp_method="none")) == 0
    rows = summary(out)
    assert rows["s1"]["Baited_reads"] == "2" and rows["s1"]["Status"] == "ok"
    assert (out / "1_extracted" / "s1.bam").exists()
    assert not list((out / "1_extracted").glob(".*.names"))


def _interrupt(monkeypatch, module, function: str, sample: str, after: bool = False):
    """Make `module.function` raise KeyboardInterrupt for `sample` (after running it, with after=True)."""
    real = getattr(module, function)

    def interrupted(name, *args, **kwargs):
        if getattr(name, "name", name) != sample:  # A sample name, or a Sample
            return real(name, *args, **kwargs)
        if after:
            real(name, *args, **kwargs)
        raise KeyboardInterrupt

    monkeypatch.setattr(module, function, interrupted)
    return lambda: monkeypatch.setattr(module, function, real)


def test_new_reads_interrupted_after_baiting_are_filtered_again(stubs, dataset, tmp_path, monkeypatch):
    import bacon.steps
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out, snp_method="none", parallel=1))
    for name in ("s1", "s2"):  # New reads: one read matches
        _reads(reads / f"{name}.fastq.gz", [("on1", 950)])
    restore = _interrupt(monkeypatch, bacon.steps, "bait_minimap2", "s2")  # After s1 was baited
    with pytest.raises(KeyboardInterrupt):
        run(settings(ref, reads, out, snp_method="none", parallel=1))
    restore()
    run(settings(ref, reads, out, snp_method="none", parallel=1))
    row = summary(out)["s1"]
    assert row["Baited_reads"] == "1" and row["Filtered_reads"] == "1"  # Not the filtered reads of the old input
    assert sum(1 for _ in read_records(out / "2_filtered" / "s1.fastq.gz")) == 1


def test_output_of_interrupted_parameter_change_is_not_reused(stubs, dataset, tmp_path, monkeypatch):
    import bacon.steps
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out, snp_method="none", parallel=1))
    restore = _interrupt(monkeypatch, bacon.steps, "filter_filtlong", "s1", after=True)
    with pytest.raises(KeyboardInterrupt):  # s1 filtered with 850, interrupted before it was recorded
        run(settings(ref, reads, out, snp_method="none", parallel=1, min_read_length=850))
    restore()
    run(settings(ref, reads, out, snp_method="none", parallel=1))  # Back to 100
    assert summary(out)["s1"]["Filtered_reads"] == "2"
    assert sum(1 for _ in read_records(out / "2_filtered" / "s1.fastq.gz")) == 2


def test_assemblies_changed_before_an_interrupted_comparison_are_compared(stubs, dataset, tmp_path, monkeypatch):
    import bacon.pipeline
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    _reads(reads / "s1.fastq.gz", [("on1", 950)])
    real = bacon.pipeline._compare
    monkeypatch.setattr(bacon.pipeline, "_compare", lambda *a, **k: (_ for _ in ()).throw(KeyboardInterrupt))
    with pytest.raises(KeyboardInterrupt):  # s1 assembled again, interrupted before the comparison
        run(settings(ref, reads, out))
    monkeypatch.setattr(bacon.pipeline, "_compare", real)
    stubs.unlink()
    run(settings(ref, reads, out))
    assert "ska" in calls(stubs)


def test_sigterm_kills_the_programs_still_running(stubs, dataset, tmp_path):
    import signal
    import subprocess
    import time
    ref, reads = dataset
    pidfile = tmp_path / "child.pid"
    flye = tmp_path / "bin" / "flye"  # A Flye that starts a program of its own and waits
    flye.write_text(f'#!/bin/sh\n[ "$1" = "--version" ] && exit 0\n'
                    f"sleep 30 & echo $! > {pidfile}.tmp; mv {pidfile}.tmp {pidfile}; wait\n")
    code = (f"import sys; from bacon.cli import main; sys.exit(main(['-r', {str(ref)!r}, '-i', {str(reads)!r}, "
            f"'-o', {str(tmp_path / 'out')!r}, '-a', 'flye', '-t', '1', '-p', '1']))")
    bacon = subprocess.Popen([sys.executable, "-c", code], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(200):
        if pidfile.exists():
            break
        time.sleep(0.05)
    child = int(pidfile.read_text())
    bacon.send_signal(signal.SIGTERM)
    assert bacon.wait(timeout=20) == 128 + signal.SIGTERM
    for _ in range(100):  # The child is gone (or a zombie of init)
        try:
            os.kill(child, 0)
        except ProcessLookupError:
            break
        if Path(f"/proc/{child}/stat").exists() and Path(f"/proc/{child}/stat").read_text().split()[2] == "Z":
            break
        time.sleep(0.05)
    else:
        os.kill(child, 9)
        raise AssertionError("the program started by Flye outlived BACoN")


def test_checkpoints_of_0_3_2_are_reused(stubs, dataset, tmp_path):
    import hashlib

    from bacon.pipeline import _fingerprint
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    md5 = hashlib.md5((out / "reference.fasta").read_bytes()).hexdigest()
    params_032 = {"reference": md5, "method": "minimap2", "kmer": None, "keep_bam": False}  # As in BACoN 0.3.2
    saved = json.loads((out / ".checkpoints" / "bait.json").read_text())
    assert saved["fingerprint"] == _fingerprint("0", params_032)


def test_a_moved_output_folder_resumes_with_its_own_files(stubs, dataset, tmp_path):
    import shutil
    ref, reads = dataset
    out, moved = tmp_path / "out", tmp_path / "moved"
    run(settings(ref, reads, out))
    shutil.move(out, moved)
    stubs.unlink()
    assert run(settings(ref, reads, moved)) == 0
    assert calls(stubs) == ["minimap2"]  # Only the failed sample is tried again
    info = json.loads((moved / "run_info.json").read_text())
    assert info["comparison"]["distances"] == str(moved / "4_compared" / "ska" / "snp_distances.tsv")


def test_output_folder_on_a_file_system_without_locks(stubs, dataset, tmp_path, monkeypatch, caplog):
    import errno
    import fcntl

    def no_locks(fd, operation):
        raise OSError(errno.ENOLCK, "No locks available")

    monkeypatch.setattr(fcntl, "flock", no_locks)
    ref, reads = dataset
    assert run(settings(ref, reads, tmp_path / "out", snp_method="none")) == 0
    assert "Could not lock" in caplog.text


def test_vcf_version_check(tmp_path):
    from bacon.pipeline import _vcf_outdated
    vcf = tmp_path / "snps.vcf"
    for source, outdated in [("BACoN 0.3.2", True), ("BACoN 0.3.3", False), ("BACoN 0.3.10", False),
                             ("BACoN 1.0", False), ("BACoN ", True), (None, True)]:
        vcf.write_text("##fileformat=VCFv4.2\n" + (f"##source={source}\n" if source else "") + "#CHROM\n")
        assert _vcf_outdated(vcf) is outdated, source


def test_a_run_after_an_interrupted_one_can_start_programs(stubs, dataset, tmp_path):
    from bacon import tools
    ref, reads = dataset
    (reads / "none.fastq.gz").unlink()  # No failed sample: the resumed run only compares
    run(settings(ref, reads, tmp_path / "out", snp_method="none"))
    tools.kill_running()  # As left by an interrupted run in this process
    assert run(settings(ref, reads, tmp_path / "out")) == 0
    assert "ska" in calls(stubs)


def test_a_moved_folder_of_0_3_3_resumes_with_its_own_files(stubs, dataset, tmp_path):
    import shutil
    ref, reads = dataset
    out, moved = tmp_path / "out", tmp_path / "moved"
    run(settings(ref, reads, out))
    for checkpoint in (out / ".checkpoints").glob("*.json"):  # As written by BACoN 0.3.3: no root
        data = json.loads(checkpoint.read_text())
        del data["root"]
        checkpoint.write_text(json.dumps(data))
    shutil.move(out, moved)
    stubs.unlink()
    assert run(settings(ref, reads, moved)) == 0
    assert calls(stubs) == ["minimap2"]
    info = json.loads((moved / "run_info.json").read_text())
    assert info["comparison"]["distances"] == str(moved / "4_compared" / "ska" / "snp_distances.tsv")


def _bacon_cli(ref, reads, out, **popen):
    """Start `bacon` (Flye assembly, one sample at a time) in a separate process."""
    import subprocess
    code = (f"import sys; from bacon.cli import main; sys.exit(main(['-r', {str(ref)!r}, '-i', {str(reads)!r}, "
            f"'-o', {str(out)!r}, '-a', 'flye', '-t', '1', '-p', '1']))")
    return subprocess.Popen([sys.executable, "-c", code], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                            **popen)


def test_nohup_runs_survive_a_closed_terminal(stubs, dataset, tmp_path):
    import signal
    import time
    ref, reads = dataset
    started = tmp_path / "started"
    flye = tmp_path / "bin" / "flye"  # A Flye that takes a while
    flye.write_text(f'#!/bin/sh\n[ "$1" = "--version" ] && exit 0\ntouch {started}\nsleep 30\n')
    bacon = _bacon_cli(ref, reads, tmp_path / "out",
                       preexec_fn=lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN))  # As nohup starts it
    for _ in range(200):
        if started.exists():
            break
        time.sleep(0.05)
    bacon.send_signal(signal.SIGHUP)  # The terminal is closed
    time.sleep(1)
    assert bacon.poll() is None  # Still running
    bacon.send_signal(signal.SIGTERM)
    assert bacon.wait(timeout=20) == 128 + signal.SIGTERM


def test_cli_in_another_thread_and_signal_handlers_restored(stubs, dataset, tmp_path):
    import signal
    import threading

    from bacon.cli import main
    ref, reads = dataset
    args = ["-r", str(ref), "-i", str(reads), "-a", "flye", "-t", "1", "--snp-method", "none"]
    codes = []
    thread = threading.Thread(target=lambda: codes.append(main([*args, "-o", str(tmp_path / "a")])))
    thread.start()
    thread.join()
    assert codes == [0]
    before = signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGHUP)
    assert main([*args, "-o", str(tmp_path / "b")]) == 0
    assert (signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGHUP)) == before


def test_errors_and_interruptions_are_in_bacon_log(stubs, dataset, tmp_path, monkeypatch):
    import bacon.steps
    ref, reads = dataset
    out = tmp_path / "out"
    with pytest.raises(BaconError, match="All samples failed"):
        run(settings(ref, reads, out, min_read_length=5000))
    assert "[ERROR] All samples failed at the filter step" in (out / "bacon.log").read_text()
    _interrupt(monkeypatch, bacon.steps, "filter_filtlong", "s1")
    with pytest.raises(KeyboardInterrupt):
        run(settings(ref, reads, out))
    assert (out / "bacon.log").read_text().rstrip().endswith("[ERROR] Interrupted")


def test_guess_root_takes_the_last_bacon_folder_in_the_path():
    from bacon.pipeline import _guess_root
    assert _guess_root({"s1": {"output": "/data/1_extracted/out/3_assembled/all_assemblies/s1.fasta"}}) == \
        "/data/1_extracted/out"
    assert _guess_root({"distances": "/x/4_compared/out/4_compared/ska/snp_distances.tsv"}) == "/x/4_compared/out"
    assert _guess_root({"s1": {"failed": "no reads"}}) is None


def _as_0_3_3(out: Path) -> None:
    """Rewrite the checkpoints of a run as BACoN 0.3.3 wrote them: no folder, no input signatures."""
    for checkpoint in (out / ".checkpoints").glob("*.json"):
        data = json.loads(checkpoint.read_text())
        data.pop("root", None)
        data["results"].pop("inputs", None)
        for res in data["results"].values():
            if isinstance(res, dict):
                res.pop("upstream", None)
        checkpoint.write_text(json.dumps(data))


def test_new_reads_interrupted_after_baiting_with_checkpoints_of_0_3_3(stubs, dataset, tmp_path, monkeypatch):
    import time

    import bacon.steps
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out, snp_method="none", parallel=1))
    _as_0_3_3(out)
    time.sleep(0.01)  # The new outputs must be newer than the old ones, whatever the file system's clock
    for name in ("s1", "s2"):
        _reads(reads / f"{name}.fastq.gz", [("on1", 950)])
    restore = _interrupt(monkeypatch, bacon.steps, "bait_minimap2", "s2")  # After s1 was baited
    with pytest.raises(KeyboardInterrupt):
        run(settings(ref, reads, out, snp_method="none", parallel=1))
    restore()
    run(settings(ref, reads, out, snp_method="none", parallel=1))
    assert summary(out)["s1"]["Filtered_reads"] == "1"
    assert sum(1 for _ in read_records(out / "2_filtered" / "s1.fastq.gz")) == 1


def test_resumed_checkpoints_of_0_3_3_get_signatures(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    _as_0_3_3(out)
    stubs.unlink()
    run(settings(ref, reads, out))
    assert calls(stubs) == ["minimap2"]  # Only the failed sample: the results of 0.3.3 are reused ...
    filt = json.loads((out / ".checkpoints" / "filter.json").read_text())
    comp = json.loads((out / ".checkpoints" / "compare.json").read_text())
    assert filt["root"] == str(out) and all("upstream" in r for r in filt["results"].values())  # ... and updated
    assert "inputs" in comp["results"]


def test_comparison_of_0_3_3_older_than_an_assembly_is_redone(stubs, dataset, tmp_path):
    import os
    import time
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    _as_0_3_3(out)
    assembly = out / "3_assembled" / "all_assemblies" / "s1.fasta"
    later = time.time() + 5
    os.utime(assembly, (later, later))  # Assembled again after the comparison (interrupted before comparing)
    stubs.unlink()
    run(settings(ref, reads, out))
    assert "ska" in calls(stubs)


def test_a_sample_with_no_read_baited_keeps_its_raw_counts(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out, snp_method="none"))
    row = summary(out)["none"]
    assert row["Raw_reads"] == "1" and row["Baited_reads"] == "0"


def test_output_folder_cannot_be_in_the_input_folder(stubs, dataset, tmp_path):
    ref, reads = dataset
    for out in (reads, reads / "bacon"):
        with pytest.raises(BaconError, match="cannot be the input folder or inside it"):
            run(settings(ref, reads, out))
    assert not (reads / "bacon").exists()


def test_guess_root_with_a_sample_named_distances():
    from bacon.pipeline import _guess_root
    assert _guess_root({"distances": {"output": "/out/1_extracted/distances.fastq.gz"}}) == "/out"


def test_memory_and_cpus_under_a_job_scheduler(tmp_path, monkeypatch):
    import bacon.pipeline as pipeline
    from bacon.cli import build_parser
    v2, v1 = tmp_path / "memory.max", tmp_path / "limit_in_bytes"
    v2.write_text("max\n")
    v1.write_text("4000000000\n")
    assert pipeline._cgroup_memory_limit((str(v2), str(v1))) == 4_000_000_000  # v2 without limit, v1 with
    assert pipeline._cgroup_memory_limit((str(tmp_path / "none"),)) is None
    monkeypatch.setattr(pipeline, "_cgroup_memory_limit", lambda: 4_000_000_000)
    assert pipeline.default_memory_gb() == 3  # 85% of 4 GB
    if hasattr(os, "sched_getaffinity"):
        monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {0, 1})
        assert pipeline.usable_cpus() == 2
        assert build_parser().parse_args(["-r", "r", "-i", "i", "-o", "o"]).threads == 2


def test_bbduk_memory_is_shared_by_the_samples_baited_in_this_run(stubs, dataset, tmp_path):
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out, baiting="bbduk", memory_gb=8, parallel=4, snp_method="none"))
    assert "-Xmx2g" in (out / "logs" / "1_bait" / "none.log").read_text()  # 4 samples at the same time
    run(settings(ref, reads, out, baiting="bbduk", memory_gb=8, parallel=4, snp_method="none"))
    assert "-Xmx8g" in (out / "logs" / "1_bait" / "none.log").read_text()  # Only the failed sample is retried


def test_run_restores_the_level_of_the_bacon_logger(stubs, dataset, tmp_path):
    import logging
    logger = logging.getLogger("bacon")
    previous = logger.level
    logger.setLevel(logging.WARNING)
    try:
        ref, reads = dataset
        run(settings(ref, reads, tmp_path / "out", snp_method="none"))
        assert logger.level == logging.WARNING
    finally:
        logger.setLevel(previous)


# ---------------------------------------------------------------------------------------------------------------
# Reference annotations
# ---------------------------------------------------------------------------------------------------------------

def _genbank_of(ref: Path, name: str = "ref", features: str = "", definition: str = "") -> str:
    """A GenBank record of the sequence of a one-record fasta, with a CDS over its first 30 bases."""
    seq = "".join(ref.read_text().split("\n")[1:])
    features = features or ("     gene            1..30\n                     /gene=\"orfA\"\n"
                            "     CDS             1..30\n                     /gene=\"orfA\"\n"
                            "                     /transl_table=11\n")
    origin = "\n".join(f"{i + 1:>9} " + " ".join(seq[j:j + 10].lower() for j in range(i, min(i + 60, len(seq)), 10))
                       for i in range(0, len(seq), 60))
    return (f"LOCUS       {name}  {len(seq)} bp    DNA     circular PLN 01-JAN-2026\n"
            + (f"DEFINITION  {definition}\n" if definition else "")
            + f"FEATURES             Location/Qualifiers\n{features}ORIGIN      \n{origin}\n//\n")


def test_genbank_reference_gives_the_same_reference_and_reruns_nothing(stubs, dataset, tmp_path):
    import hashlib
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    before = (out / "reference.fasta").read_bytes()
    checkpoints = {p.name: json.loads(p.read_text())["fingerprint"] for p in (out / ".checkpoints").glob("*.json")}
    gb = tmp_path / "ref.gb"
    gb.write_text(_genbank_of(ref))
    stubs.unlink()
    assert run(settings(gb, reads, out)) == 0
    assert calls(stubs) == ["minimap2"]  # Only the sample that failed before: no step ran again
    assert (out / "reference.fasta").read_bytes() == before
    assert {p.name: json.loads(p.read_text())["fingerprint"] for p in (out / ".checkpoints").glob("*.json")} == \
        checkpoints
    info = json.loads((out / "run_info.json").read_text())
    assert info["reference"]["md5"] == hashlib.md5(gb.read_bytes()).hexdigest()  # The file given
    assert info["annotation"]["format"] == "genbank" and info["annotation"]["copy"] == "annotation.gb"
    assert info["annotation"]["genes"] == 1 and info["settings"]["annotation"] is None
    assert (out / "annotation.gb").read_text() == gb.read_text()
    page = (out / "report.html").read_text()
    # The VCF's SNP at 5 (A>C) is the middle base of codon 2 (TAC): Y2S
    assert 'class="gene-cds"' in page and ">Y2S<" in page and ">missense<" in page and "orfA" in page
    assert "translation table 11" in page
    # Back to the fasta: still nothing to rerun, and the annotation copy goes
    stubs.unlink()
    assert run(settings(ref, reads, out)) == 0
    assert calls(stubs) == ["minimap2"] and not (out / "annotation.gb").exists()
    assert json.loads((out / "run_info.json").read_text())["annotation"] is None


def test_annotation_option_is_not_in_the_checkpoints(stubs, dataset, tmp_path):
    import gzip
    import shutil
    ref, reads = dataset
    out = tmp_path / "out"
    run(settings(ref, reads, out))
    gb = tmp_path / "ann.gbk.gz"
    with gzip.open(gb, "wt") as fh:
        fh.write(_genbank_of(ref))
    stubs.unlink()
    assert run(settings(ref, reads, out, annotation=gb)) == 0
    assert calls(stubs) == ["minimap2"]
    copy = out / "annotation.gb"
    assert copy.read_text().startswith("LOCUS")  # Decompressed
    info = json.loads((out / "run_info.json").read_text())
    assert info["settings"]["annotation"] == str(gb) and info["annotation"]["file"] == str(gb)
    assert ">Y2S<" in (out / "report.html").read_text()
    # A moved folder: the report is rebuilt from the copy
    moved = tmp_path / "moved"
    shutil.move(out, moved)
    from bacon.report import write_report
    assert ">Y2S<" in write_report(moved).read_text()
    # A GFF3 annotation replaces the GenBank copy
    gff = tmp_path / "ann.gff3"
    gff.write_text("##gff-version 3\nref\tt\tgene\t1\t30\t.\t+\t.\tID=g1;gene=orfB\n"
                   "ref\tt\tCDS\t1\t30\t.\t+\t0\tID=c1;Parent=g1\n")
    stubs.unlink()
    assert run(settings(ref, reads, moved, annotation=gff)) == 0
    assert calls(stubs) == ["minimap2"]
    assert (moved / "annotation.gff3").exists() and not (moved / "annotation.gb").exists()
    assert "orfB" in (moved / "report.html").read_text()


def test_annotation_errors_and_warnings(stubs, dataset, tmp_path, caplog):
    ref, reads = dataset
    out = tmp_path / "out"
    with pytest.raises(BaconError, match="Annotation file not found"):
        run(settings(ref, reads, out, annotation=tmp_path / "missing.gb", snp_method="none"))
    with pytest.raises(BaconError, match="not a GenBank"):
        run(settings(ref, reads, out, annotation=ref, snp_method="none"))
    assert not (out / "1_extracted").exists()  # Failed before any step
    other = tmp_path / "other.gb"
    other.write_text(_genbank_of(ref, name="other").replace("1000 bp", "999 bp"))
    assert run(settings(ref, reads, out, annotation=other, snp_method="none")) == 0
    assert "match no reference sequence by name" in caplog.text
    beyond = tmp_path / "beyond.gb"
    beyond.write_text(_genbank_of(ref, features="     gene            900..1200\n                     /gene=\"x\"\n"))
    caplog.clear()
    assert run(settings(ref, reads, out, annotation=beyond, snp_method="none")) == 0
    assert "1 feature(s) beyond the end" in caplog.text


def test_gff3_cannot_be_the_reference(stubs, dataset, tmp_path):
    ref, reads = dataset
    gff = tmp_path / "ref.gff3"
    gff.write_text("##gff-version 3\nref\tt\tgene\t1\t30\t.\t+\t.\tID=g1\n")
    with pytest.raises(BaconError, match="not GFF3"):
        run(settings(gff, reads, tmp_path / "out", snp_method="none"))
    empty = tmp_path / "noseq.gb"
    empty.write_text(_genbank_of(ref).split("ORIGIN")[0] + "//\n")
    with pytest.raises(BaconError, match="no sequence"):
        run(settings(empty, reads, tmp_path / "out", snp_method="none"))


def test_genbank_reference_with_a_definition_matches_ncbi_fasta_headers(stubs, dataset, tmp_path):
    ref, reads = dataset
    gb = tmp_path / "NC_1.gb"
    gb.write_text(_genbank_of(ref, name="NC_1", definition="Some plant chloroplast, complete genome.")
                  .replace("LOCUS       NC_1 ", "LOCUS       NC_000001 ").replace("FEATURES", "VERSION     NC_1.2\nFEATURES"))
    out = tmp_path / "out"
    assert run(settings(gb, reads, out, snp_method="none")) == 0
    assert (out / "reference.fasta").read_text().startswith(">NC_1.2 Some plant chloroplast, complete genome\nACGT")
