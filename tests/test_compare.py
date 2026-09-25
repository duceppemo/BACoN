from bacon.compare import clean_alignment, snp_distances, write_distances
from bacon.seqio import Record, read_records


def test_clean_alignment_and_rename(tmp_path):
    src = tmp_path / "a.fasta"
    src.write_text(">x.fasta\nacgtRN-\n>ref.fasta.ref\nACGTAAA\n")
    recs = clean_alignment(src, tmp_path / "b.fasta", {"x.fasta": "x", "ref.fasta.ref": "Reference"})
    assert [(r.name, r.seq) for r in recs] == [("x", "ACGTNN-"), ("Reference", "ACGTAAA")]
    assert [r.name for r in read_records(tmp_path / "b.fasta")] == ["x", "Reference"]


def test_snp_distances_ignore_gaps_and_n_and_sort_reference_first(tmp_path):
    recs = [Record("b", "ACGTN-A"), Record("Reference", "ACGTAAA"), Record("a", "TCGAAAC")]
    names, m = snp_distances(recs)
    assert names == ["Reference", "a", "b"]
    assert m == [[0, 3, 0], [3, 0, 3], [0, 3, 0]]
    out = tmp_path / "d.tsv"
    write_distances(out, names, m)
    assert out.read_text().splitlines()[0] == "snp-dists\tReference\ta\tb"


def test_wrap_circular():
    from bacon.compare import wrap_circular
    circ = Record("s_1 circular=true", "ACGTACGTAC")
    assert wrap_circular(circ, 4).seq == "ACGTACGTACACG"
    assert wrap_circular(Record("s_2", "ACGTACGTAC"), 4).seq == "ACGTACGTAC"
    assert wrap_circular(Record("ref", "ACGTACGTAC"), 4, force=True).seq == "ACGTACGTACACG"
    assert wrap_circular(Record("tiny circular=true", "ACG"), 4).seq == "ACG"
