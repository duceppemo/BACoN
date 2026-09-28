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


def _vcf(tmp_path, body, samples=("Reference", "a", "b"), header=True):
    raw = tmp_path / "raw.vcf"
    head = "##fileformat=VCFv4.4\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
    raw.write_text((head if header else "") + body)
    return raw


def _clean(tmp_path, raw):
    from bacon.compare import clean_vcf
    ref = tmp_path / "ref.fasta"
    ref.write_text(">chr\n" + "A" * 100 + "\n")
    out = tmp_path / "out.vcf"
    n = clean_vcf(raw, out, rename={}, reference=ref, source="t")
    return n, [line.split("\t") for line in out.read_text().splitlines() if not line.startswith("#")]


def test_clean_vcf_n_is_not_an_allele(tmp_path):
    raw = _vcf(tmp_path, "chr\t5\t.\tA\tN\t.\t.\t.\tGT\t0\t1\t0\n"          # N only: dropped
                         "chr\t6\t.\tA\tG,T,N\t40\tN\tNA\tGT\t0\t3\t2\n"   # N removed, 3 -> ., 2 stays 2
                         "chr\t7\t.\tA\tN,C\t40\tN;IND\tNA\tGT\t0\t1\t2\n"  # C renumbered to 1
                         "chr\t8\t.\tA\tC\t.\t.\t.\tGT\t1\t1\t1\n"          # reference not 0: dropped
                         "\n")
    n, rows = _clean(tmp_path, raw)
    assert n == 2
    assert rows[0][1:] == ["6", ".", "A", "G,T", "40", "PASS", ".", "GT", ".", "2"]
    assert rows[1][1:] == ["7", ".", "A", "C", "40", "IND", ".", "GT", ".", "1"]


def test_clean_vcf_folds_positions_of_an_extended_circular_reference(tmp_path):
    raw = _vcf(tmp_path, "chr\t30\t.\tA\tC\t.\t.\t.\tGT\t0\t1\t0\n"
                         "chr\t103\t.\tA\tG\t.\t.\t.\tGT\t0\t0\t1\n"   # beyond the 100 bp contig: position 3
                         "chr\t105\t.\tA\tT\t.\t.\t.\tGT\t0\t1\t0\n"
                         "chr\t5\t.\tA\tT\t.\t.\t.\tGT\t0\t0\t1\n")    # duplicate of 105 folded: dropped
    n, rows = _clean(tmp_path, raw)
    assert n == 3 and [r[1] for r in rows] == ["3", "5", "30"]
    assert rows[1][9:] == ["1", "0"]  # The first record at a position is kept


def test_clean_vcf_rejects_malformed_input(tmp_path):
    import pytest

    from bacon import BaconError
    for body, header, samples, message in [
        ("chr\t5\t.\tA\tC\n", True, ("Reference", "a"), "truncated"),
        ("chr\t5\t.\tA\tC\t.\t.\t.\tGT\t0\t1\n", False, ("Reference", "a"), "before the #CHROM"),
        ("", False, ("Reference", "a"), "no #CHROM"),
        ("", True, ("a", "b"), "one reference column"),
    ]:
        raw = _vcf(tmp_path, body, samples=samples, header=header)
        with pytest.raises(BaconError, match=message):
            _clean(tmp_path, raw)
    assert not (tmp_path / "out.vcf.tmp").exists()


def test_ska_vcf_has_snps_next_to_the_ends_of_circular_genomes(tmp_path):
    """With the real SKA2: SNPs within 15 bp of the start of a circular reference are in the VCF as in the
    alignment (the reference is extended for both)."""
    import random
    import shutil

    import pytest
    if not shutil.which("ska"):
        pytest.skip("ska not installed")
    from bacon.compare import run_ska, write_vcf
    rng = random.Random(1)
    ref = "".join(rng.choice("ACGT") for _ in range(5000))
    (tmp_path / "ref.fasta").write_text(f">chr\n{ref}\n")
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    assemblies = {}
    for name, pos in {"a": 5, "b": 4980, "c": 2000, "d": None}.items():
        g = list(ref)
        if pos is not None:
            g[pos] = swap[g[pos]]
        g = "".join(g)
        path = tmp_path / f"{name}.fasta"
        path.write_text(f">{name}_1 circular=true\n{g[1000:] + g[:1000]}\n")  # Rotated, like an assembly
        assemblies[name] = path
    (tmp_path / "logs").mkdir()
    run_ska(tmp_path / "ref.fasta", assemblies, tmp_path / "out", tmp_path / "logs", threads=1, min_freq=1.0)
    vcf, n = write_vcf("ska", tmp_path / "ref.fasta", tmp_path / "out", tmp_path / "logs", threads=1,
                       assemblies=assemblies, source="t")
    rows = [line.split("\t") for line in vcf.read_text().splitlines() if not line.startswith("#")]
    assert [(r[1], r[9:]) for r in rows] == [("6", ["1", "0", "0", "0"]), ("2001", ["0", "0", "1", "0"]),
                                            ("4981", ["0", "1", "0", "0"])]
    assert not (tmp_path / "out" / "snps.raw.vcf").exists()
