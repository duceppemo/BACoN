import pytest

from bacon import BaconError
from bacon.seqio import (
    check_reference,
    extract_reads,
    n50,
    read_records,
    read_stats,
    sniff_format,
    split_extension,
)


@pytest.mark.parametrize("name, expected", [
    ("s1.fastq.gz", ("s1", "fastq")), ("s1.fq", ("s1", "fastq")), ("a.b.fasta", ("a.b", "fasta")),
    ("x.FNA.GZ", ("x", "fasta")), ("notes.txt", None), (".fastq", None), ("s.fastq.bz2", None),
])
def test_split_extension(name, expected):
    assert split_extension(name) == expected


def test_read_fastq_gz_and_plain(fastq):
    for gz in (True, False):
        path = fastq(f"r{gz}.fastq" + (".gz" if gz else ""), [("r1", "ACGT"), ("r2", "AAAAAA")], gz)
        recs = list(read_records(path))
        assert [(r.name, r.seq) for r in recs] == [("r1", "ACGT"), ("r2", "AAAAAA")]
        assert recs[0].header == "r1 extra"


def test_gzip_detected_by_content_not_name(tmp_path, fastq):
    path = fastq("mislabelled.fastq", [("r1", "ACGT")], gz=True)  # gzipped, no .gz
    assert [r.seq for r in read_records(path)] == ["ACGT"]


def test_multiline_fasta(tmp_path):
    p = tmp_path / "m.fasta"
    p.write_text(">a desc\nACG\nTT\n\n>b\nGG\n")
    assert [(r.name, r.seq) for r in read_records(p)] == [("a", "ACGTT"), ("b", "GG")]


def test_truncated_fastq_is_an_error(tmp_path):
    p = tmp_path / "t.fastq"
    p.write_text("@r1\nACGT\n+\nII\n")
    with pytest.raises(BaconError, match="malformed or truncated"):
        list(read_records(p))


def test_sniff_format(tmp_path):
    (tmp_path / "e.fastq").write_text("\n\n")
    assert sniff_format(tmp_path / "e.fastq") is None
    (tmp_path / "x.fasta").write_text("hello\n")
    with pytest.raises(BaconError, match="neither fasta nor fastq"):
        sniff_format(tmp_path / "x.fasta")


def test_n50():
    assert n50([]) == 0
    assert n50([10]) == 10
    assert n50([2, 3, 4, 10]) == 10  # 10 >= 19 / 2
    assert n50([5, 5, 5, 5]) == 5
    assert n50([1, 2, 3, 4, 5]) == 4  # 5 + 4 = 9 >= 7.5


def test_read_stats(fastq):
    stats = read_stats([fastq("a.fastq.gz", [("r1", "ACGT"), ("r2", "A" * 10)])])
    assert (stats.reads, stats.bases, stats.n50) == (2, 14, 10)


def test_extract_reads(fastq, tmp_path):
    src1 = fastq("a.fastq.gz", [("r1", "ACGT"), ("r2", "GGGG")])
    src2 = fastq("b.fastq", [("r3", "TTTTTT")], gz=False)
    out = tmp_path / "out.fastq.gz"
    total, kept = extract_reads([src1, src2], {"r2", "r3", "absent"}, out)
    assert (total.reads, total.bases, kept.reads, kept.bases) == (3, 14, 2, 10)
    assert [r.name for r in read_records(out)] == ["r2", "r3"]
    assert next(read_records(out)).qual == "IIII"


def test_check_reference(fasta):
    assert check_reference(fasta("ok.fasta", [("a", "ACGT"), ("b", "GG")])) == [("a", 4), ("b", 2)]
    with pytest.raises(BaconError, match="duplicate"):
        check_reference(fasta("dup.fasta", [("a", "ACGT"), ("a", "GG")]))
    with pytest.raises(BaconError, match="no sequence"):
        check_reference(fasta("empty.fasta", [("a", "")]))


def test_check_reference_rejects_fastq_and_missing(fastq, tmp_path):
    with pytest.raises(BaconError, match="must be a fasta"):
        check_reference(fastq("r.fastq.gz", [("r1", "ACGT")]))
    with pytest.raises(BaconError, match="not found"):
        check_reference(tmp_path / "missing.fasta")
