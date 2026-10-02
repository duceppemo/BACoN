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
