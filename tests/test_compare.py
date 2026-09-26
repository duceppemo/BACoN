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


def test_clean_vcf_renames_drops_reference_and_adds_contigs(tmp_path):
    from bacon.compare import clean_vcf
    ref = tmp_path / "ref.fasta"
    ref.write_text(">chr1\nACGTACGT\n>chr2\nAC\n")
    raw = tmp_path / "raw.vcf"
    raw.write_text("##fileformat=VCFv4.4\n##contig=<ID=chr1>\n##source=harvest\n"
                   '##FILTER=<ID=IND,Description="indel">\n'
                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tref.fasta.ref\ta.fasta\tb.fasta\n"
                   "chr1\t3\t.\tG\tT\t40\tPASS\tNA\tGT\t0\t1\t0\n"
                   "chr1\t5\t.\tA\tC\t40\tIND\tNA\tGT\t0\t0\t1\n"
                   "chr1\t6\t.\tC\t.\t.\t.\t.\tGT\t0\t.\t0\n")  # No alternate allele: dropped
    out = tmp_path / "out.vcf"
    n = clean_vcf(raw, out, rename={"ref.fasta.ref": "Reference", "a.fasta": "a", "b.fasta": "b"}, reference=ref,
                  source="BACoN test")
    lines = out.read_text().splitlines()
    assert n == 2
    assert lines[:4] == ["##fileformat=VCFv4.2", "##source=BACoN test", "##contig=<ID=chr1,length=8>",
                         "##contig=<ID=chr2,length=2>"]
    assert '##FILTER=<ID=IND,Description="indel">' in lines and "##source=harvest" not in lines
    assert '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' in lines
    assert lines[-3].split("\t")[9:] == ["a", "b"]
    assert lines[-2].split("\t")[9:] == ["1", "0"] and lines[-1].split("\t")[6] == "IND"
