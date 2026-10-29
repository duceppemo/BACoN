import shutil

import pytest

from bacon.steps import _fill_from, _flye_info, count_circular


def test_count_circular_headers(tmp_path):
    p = tmp_path / "a.fasta"
    p.write_text(">u1ctg_len-70000_circular-yes_depth-30 mult=1\nA\n"   # myloasm
                 ">u2ctg_len-9000_circular-no_depth-30\nA\n"
                 ">Utg30 LN:i:69995 RC:i:15 XO:i:1\nA\n"                  # Raven
                 ">Utg31 LN:i:100 XO:i:0\nA\n"
                 ">contig_1 circular=true\nA\n>plain\nA\n")
    assert count_circular(p) == 3


def test_flye_info(tmp_path):
    p = tmp_path / "assembly_info.txt"
    p.write_text("#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
                 "contig_1\t1000\t30\tY\tN\t1\t*\t1\ncontig_2\t500\t10\tN\tN\t1\t*\t2\n")
    info = _flye_info(p)
    assert info["contig_1"]["circ."] == "Y" and info["contig_2"]["length"] == "500"
    assert _flye_info(tmp_path / "missing.txt") == {}


def test_fill_from_reference_only_at_the_ends():
    assert _fill_from("NNACGTNNTN", "GGACGTAAAA") == "GGACGTNNTA"
    assert _fill_from("ACGT", "") == "ACGT"
    assert _fill_from("NNNN", "ACGT") == "NNNN"


def test_publish_marks_circular_contigs(tmp_path):
    from bacon.seqio import read_records
    from bacon.steps import AssemblyDirs, _publish_assembly
    raw = tmp_path / "raw.fasta"
    raw.write_text(">contig_1\nACGT\n>u2ctg_circular-yes_depth-3\nAC\n>contig_3\nGG\n")
    dirs = AssemblyDirs(tmp_path / "asm")
    _publish_assembly("s1", raw, dirs, rename=False, circular={"contig_3"})
    headers = [r.header for r in read_records(dirs.assemblies / "s1.fasta")]
    assert headers == ["s1_contig_1", "s1_u2ctg_circular-yes_depth-3 circular=true", "s1_contig_3 circular=true"]


def test_filter_by_length_keeps_longest_up_to_target(tmp_path):
    from bacon.seqio import read_records
    from bacon.steps import filter_by_length
    src = tmp_path / "r.fasta"
    src.write_text("".join(f">r{i}\n{'A' * n}\n" for i, n in enumerate([100, 900, 50, 700, 700, 300])))
    out = tmp_path / "o.fasta.gz"
    filter_by_length(src, out, min_length=80, target_bases=1_800)
    assert [r.name for r in read_records(out)] == ["r1", "r3", "r4"]  # 900 + 700 + 700 >= 1,800
    filter_by_length(src, out, min_length=80, target_bases=1_500)
    assert [r.name for r in read_records(out)] == ["r1", "r3"]  # 900 + 700 reach 1,500: one of the two 700s
    filter_by_length(src, out, min_length=80, target_bases=10_000)
    assert [r.name for r in read_records(out)] == ["r0", "r1", "r3", "r4", "r5"]  # All but the 50 bp read


def test_filter_by_length_keeps_nothing_when_no_read_is_long_enough(tmp_path):
    from bacon.seqio import read_records
    from bacon.steps import filter_by_length
    src = tmp_path / "r.fasta"
    src.write_text(">r1\nACGT\n>r2\nACGTACGT\n")
    out = tmp_path / "o.fasta.gz"
    filter_by_length(src, out, min_length=100, target_bases=1000)
    assert list(read_records(out)) == []


def test_pipe_failure_names_the_program_that_failed(tmp_path):
    import pytest

    from bacon.tools import ToolError, run
    with pytest.raises(ToolError, match="^false failed with exit code 1"):
        run([["yes"], ["false"]], tmp_path / "log")
    with pytest.raises(ToolError, match="^sh failed with exit code 7"):
        run([["sh", "-c", "exit 7"], ["cat"]], tmp_path / "log")


def test_published_assemblies_are_upper_case_acgtn(tmp_path):
    from bacon.seqio import read_records
    from bacon.steps import AssemblyDirs, _publish_assembly
    raw = tmp_path / "raw.fasta"
    raw.write_text(">c1\nacgtRYKMn\n")
    dirs = AssemblyDirs(tmp_path / "asm")
    _publish_assembly("s1", raw, dirs, rename=False)
    assert next(read_records(dirs.assemblies / "s1.fasta")).seq == "ACGTNNNNN"


def test_no_program_starts_once_bacon_is_stopping(tmp_path):
    import pytest

    from bacon.tools import ToolError, allow_programs, kill_running, run
    kill_running()
    try:
        with pytest.raises(ToolError, match="not started: BACoN is stopping"):
            run(["true"], tmp_path / "log")
    finally:
        allow_programs()
    run(["true"], tmp_path / "log")


def test_a_second_interruption_during_the_kill_leaves_no_program_running(tmp_path, monkeypatch):
    import signal
    import subprocess
    import threading

    import pytest

    from bacon import tools
    real_kill, seen = tools._kill, []

    def kill(proc):  # A second Ctrl-C arrives once the first program is killed, before the second one is
        seen.append(proc)
        if len(seen) == 2:
            raise KeyboardInterrupt
        real_kill(proc)

    monkeypatch.setattr(tools, "_kill", kill)
    timer = threading.Timer(0.5, signal.pthread_kill, (threading.main_thread().ident, signal.SIGINT))
    timer.start()  # The first Ctrl-C
    try:
        with pytest.raises(KeyboardInterrupt):
            tools.run([["sleep", "30"], ["sleep", "30"]], tmp_path / "log")
    finally:
        timer.cancel()
        monkeypatch.setattr(tools, "_kill", real_kill)
    assert len(seen) == 2
    tools.kill_running()  # As BACoN does at the end of an interrupted run
    for proc in seen:
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            real_kill(proc)
            raise AssertionError(f"{proc.args} still running") from None


def test_templated_assembly_notes_reference_sequences_without_reads(tmp_path, monkeypatch):
    from pathlib import Path

    import bacon.steps as steps

    def fake_run(cmd, log_file, **kwargs):  # The consensus has only the reference sequence with reads
        if "consensus" in cmd:
            Path(cmd[cmd.index("-o") + 1]).write_text(">chr\nACGTACGTAC\n")

    monkeypatch.setattr(steps, "run", fake_run)
    reference = tmp_path / "ref.fasta"
    reference.write_text(">chr\nACGTACGTAC\n>plasmid\nGGGGCCCC\n")
    res = steps.assemble_samtools("s1", tmp_path / "reads.fq", reference, steps.AssemblyDirs(tmp_path / "asm"),
                                  tmp_path / "logs", fill_gaps=False, threads=1, reference_length=18)
    assert "no reads on reference sequence(s) plasmid: not in the assembly" in res.notes


def test_only_a_sample_folder_in_the_step_folder_is_deleted(tmp_path):
    import pytest

    from bacon.steps import AssemblyDirs, SampleFailed
    root = tmp_path / "out" / "3_assembled"
    (root / "s1").mkdir(parents=True)
    (root / "s1" / "old.fasta").write_text(">x\nA\n")
    (tmp_path / "out" / "reads.fastq").write_text("@r\nA\n+\nI\n")
    dirs = AssemblyDirs(root)
    assert dirs.clear_work("s1") == root / "s1" and not (root / "s1").exists()
    for name in ("..", ".", "s1/.."):
        with pytest.raises(SampleFailed, match="refusing to delete"):
            dirs.clear_work(name)
    assert (tmp_path / "out" / "reads.fastq").is_file() and root.is_dir()


def test_logged_command_line_can_be_run_again(tmp_path):
    """The working folder of a program is in its logged command line."""
    from bacon.tools import run
    work = tmp_path / "my work"
    work.mkdir()
    run(["sh", "-c", "pwd"], tmp_path / "log", cwd=work)
    run(["echo", "a b"], tmp_path / "log", stdout=tmp_path / "out put.txt")
    lines = (tmp_path / "log").read_text().splitlines()
    assert lines[:2] == [f"$ (cd '{work}' && sh -c pwd)", str(work)]
    assert lines[2] == f"$ echo 'a b' > '{tmp_path / 'out put.txt'}'"


@pytest.mark.skipif(not shutil.which("bbduk.sh"), reason="BBDuk is not installed")
@pytest.mark.parametrize("names", [["v.FASTQ.GZ"], ["y.FQ.GZ"], ["z.fastq"], ["x.fastq.gz"], ["w.fq.gz"],
                                   ["a.FQ.GZ", "b.fastq"]])
def test_bbduk_reads_files_whatever_their_extension_says(tmp_path, names):
    """BBTools picks the format and compression from the extension (case-sensitively); BACoN from the content."""
    import gzip
    import random

    from bacon.samples import Sample
    from bacon.steps import bait_bbduk
    rng = random.Random(1)
    genome = "".join(rng.choice("ACGT") for _ in range(3000))
    ref = tmp_path / "ref.fasta"
    ref.write_text(f">r\n{genome}\n")
    files = []
    for i, name in enumerate(names):
        reads = "".join(f"@r{i}_{j}\n{genome[j * 100:j * 100 + 300]}\n+\n{'I' * 300}\n" for j in range(5))
        path = tmp_path / name
        # Gzipped or not, against what the name says, for the two lower-case cases
        compressed = name.upper().endswith(".GZ") != (name in ("z.fastq", "x.fastq.gz"))
        path.write_bytes(gzip.compress(reads.encode()) if compressed else reads.encode())
        files.append(path)
    result = bait_bbduk(Sample("s", files, "fastq"), ref, tmp_path / "out", tmp_path, threads=1, kmer=31,
                        hdist=0, memory_gb=1)
    assert result.stats["Baited_reads"] == 5 * len(names)
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == ["s.fastq.gz"]  # No link left behind


@pytest.mark.skipif(not shutil.which("minimap2"), reason="minimap2 is not installed")
def test_off_reference_share(tmp_path):
    import random

    from bacon.steps import off_reference
    rng = random.Random(1)
    ref = "".join(rng.choice("ACGT") for _ in range(20000))
    other = "".join(rng.choice("ACGT") for _ in range(5000))
    (tmp_path / "ref.fa").write_text(f">ref\n{ref}\n")
    (tmp_path / "asm.fa").write_text(f">c1 circular=true\n{ref[:15000]}\n>c2\n{other}\n"
                                     f">c3\n{ref[15000:]}{other[:2000]}\n")
    pct = off_reference(tmp_path / "asm.fa", tmp_path / "ref.fa", tmp_path / "log", 1)
    assert 25 < pct < 27  # 7000 of 27000 bases
    (tmp_path / "same.fa").write_text(f">c1\n{ref}\n")
    assert off_reference(tmp_path / "same.fa", tmp_path / "ref.fa", tmp_path / "log", 1) < 1


def test_spans_from_paf_sums_split_alignments(tmp_path):
    from bacon.steps import _spans_from_paf
    paf = tmp_path / "a.paf"
    paf.write_text("r1\t1000\t0\t400\t+\tref\t9000\t0\t400\t380\t400\t60\n"
                   "r1\t1000\t600\t1000\t+\tref\t9000\t8600\t9000\t390\t400\t60\n"   # Supplementary piece
                   "r2\t20000\t100\t300\t+\tref\t9000\t50\t250\t180\t200\t60\n")
    spans = _spans_from_paf(paf)
    assert spans == {"r1": (1000, 800), "r2": (20000, 200)}


def test_spans_from_sam_names_cigars(tmp_path):
    from bacon.steps import _spans_from_sam_names_cigars
    lines = tmp_path / "names"
    # r1: primary 400M600S, supplementary 600H400M (hard clips keep the full length); r2: mostly clipped
    lines.write_text("r1\t400M600S\nr1\t600H400M\nr2\t100S200M19700S\nr3\t*\n")
    spans = _spans_from_sam_names_cigars(lines)
    assert spans == {"r1": (1000, 800), "r2": (20000, 200)}


def test_mostly_on_reference_fraction_of_read_or_reference():
    from bacon.steps import mostly_on_reference
    spans = {"mito": (1000, 800),        # 80% of the read aligned
             "numt": (20000, 200),       # A 200 bp homologous insert in a long nuclear read
             "locus": (20000, 8500)}     # Covers a 9 kb reference whole: passes, however long the read
    assert mostly_on_reference(spans, 9000, 0.5) == {"mito", "locus"}
    assert mostly_on_reference(spans, 9000, 0.0) == {"mito", "numt", "locus"}


@pytest.mark.skipif(not shutil.which("minimap2"), reason="minimap2 is not installed")
def test_bait_minimap2_drops_numt_like_reads(tmp_path):
    import random

    from bacon.samples import Sample
    from bacon.steps import bait_minimap2
    rng = random.Random(7)
    genome = "".join(rng.choice("ACGT") for _ in range(3000))
    other = "".join(rng.choice("ACGT") for _ in range(3000))
    ref = tmp_path / "ref.fasta"
    ref.write_text(f">r\n{genome}\n")
    reads = tmp_path / "reads.fastq"
    on_target = genome[:1500]
    numt_like = other[:1200] + genome[500:1000] + other[1200:]  # 500 of 4,700 bp on the reference
    reads.write_text(f"@mito\n{on_target}\n+\n{'I' * len(on_target)}\n"
                     f"@numt\n{numt_like}\n+\n{'I' * len(numt_like)}\n")
    sample = Sample("s", [reads], "fastq")
    strict = bait_minimap2(sample, ref, tmp_path / "strict", tmp_path, threads=1, keep_bam=False,
                           reference_length=3000, min_aligned_frac=0.5)
    assert strict.stats["Baited_reads"] == 1 and strict.stats["Bait_hit_reads"] == 2
    loose = bait_minimap2(sample, ref, tmp_path / "loose", tmp_path, threads=1, keep_bam=False,
                          reference_length=3000, min_aligned_frac=0.0)
    assert loose.stats["Baited_reads"] == 2


def test_assemble_flye_meta_flag(tmp_path, monkeypatch):
    from pathlib import Path

    import bacon.steps as steps_module
    from bacon.steps import AssemblyDirs, assemble_flye
    commands = []

    def fake_run(cmd, log_file, **kwargs):
        commands.append(cmd)
        if cmd[0] == "flye":
            out = Path(cmd[cmd.index("--out-dir") + 1])
            out.mkdir(parents=True, exist_ok=True)
            (out / "assembly.fasta").write_text(">contig_1\nACGT\n")
            (out / "assembly_info.txt").write_text(
                "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
                "contig_1\t4\t30\tY\tN\t1\t*\t1\n")

    monkeypatch.setattr(steps_module, "run", fake_run)
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    for meta, expected in ((True, True), (False, False)):
        dirs = AssemblyDirs(tmp_path / f"asm_{meta}")
        result = assemble_flye("s", reads, dirs, tmp_path, genome_size=4, read_type="nano-hq",
                               min_overlap=None, iterations=1, threads=1, reference_length=4, meta=meta)
        assert ("--meta" in commands[-1]) is expected
        assert result.stats["Circular_contigs"] == 1
