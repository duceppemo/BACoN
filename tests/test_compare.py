from pathlib import Path

import pytest

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
                         "chr\t9\t.\tN\tA\t40\tPASS\tNA\tGT\t0\t1\t1\n"     # reference base N (HarvestTools): dropped
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


def test_genotype_renumbering_handles_diploid_phased_and_extra_fields():
    from bacon.compare import _genotype_fix
    remap = {"1": ".", "2": "1"}
    assert _genotype_fix("2", remap) == "1"
    assert _genotype_fix("0/2", remap) == "0/1"
    assert _genotype_fix("1|2:35", remap) == ".|1:35"
    assert _genotype_fix(".", remap) == "."


def test_write_ska_reference_extends_when_most_assemblies_are_circular(tmp_path):
    from bacon.compare import write_ska_reference
    ref = tmp_path / "ref.fasta"
    ref.write_text(">chr\n" + "ACGT" * 20 + "\n")
    asm = {}
    for i, circular in enumerate([True, True, False]):
        asm[str(i)] = tmp_path / f"{i}.fasta"
        asm[str(i)].write_text(f">c{i}{' circular=true' if circular else ''}\nACGT\n")
    assert write_ska_reference(ref, asm, tmp_path / "out.fasta", kmer=5)
    assert next(read_records(tmp_path / "out.fasta")).seq == "ACGT" * 20 + "ACGT"
    assert not write_ska_reference(ref, {"0": asm["2"]}, tmp_path / "out.fasta", kmer=5)


def test_parsnp_warns_when_the_core_is_a_fraction_of_the_reference(tmp_path, monkeypatch, caplog):
    import bacon.compare as compare
    ref = tmp_path / "reference.fasta"
    ref.write_text(">chr\n" + "A" * 1000 + "\n")
    asm = {"a": tmp_path / "a.fasta"}
    asm["a"].write_text(">a\n" + "A" * 300 + "\n")

    def fake_run(cmd, log, **kwargs):
        if cmd[0] == "parsnp":
            assert not any(" " in arg for arg in cmd)  # Parsnp splits its command lines on spaces
            out = kwargs["cwd"] / cmd[cmd.index("-o") + 1]
            out.mkdir(parents=True)
            assert [(kwargs["cwd"] / f).read_text() for f in cmd[cmd.index("-d") + 1:cmd.index("-o")]] == [
                asm["a"].read_text()]
            (out / "parsnp.xmfa").write_text("x")
            (out / "parsnp.ggr").write_text("x")
        elif "-M" in cmd:
            Path(cmd[cmd.index("-M") + 1]).write_text(">reference.fasta.ref\n" + "A" * 300 + "\n>a.fasta\n" + "A" * 300 + "\n")
        elif "-S" in cmd:
            Path(cmd[cmd.index("-S") + 1]).write_text(">reference.fasta.ref\nAA\n>a.fasta\nCA\n")  # SNP, constant

    monkeypatch.setattr(compare, "run", fake_run)
    (tmp_path / "logs").mkdir()
    core, snps = compare.run_parsnp(ref, asm, tmp_path / "out", tmp_path / "logs", threads=1)
    assert "core genome is 300 bp, 30% of the reference" in caplog.text
    assert [r.seq for r in read_records(snps)] == ["A", "C"]  # SNP sites only
    assert {len(r.seq) for r in read_records(core)} == {300}  # The core genome keeps every column


def test_circular_rule_counts_samples_not_added_genomes(tmp_path):
    from bacon.compare import circular_genomes
    asm = {}
    for name, circular in [("a", True), ("b", True), ("pub1", False), ("pub2", False), ("pub3", False)]:
        asm[name] = tmp_path / f"{name}.fasta"
        asm[name].write_text(f">{name}{' circular=true' if circular else ''}\nACGT\n")
    assert circular_genomes(asm, added={"pub1", "pub2", "pub3"})  # 2 of 2 samples
    assert not circular_genomes(asm)  # 2 of 5 genomes
    assert not circular_genomes({k: asm[k] for k in ("pub1",)}, added={"pub1"})  # No sample at all


def test_added_genomes_are_extended_like_the_reference(tmp_path, monkeypatch):
    import bacon.compare as compare
    seen = {}

    def fake_run(cmd, log, **kwargs):
        if cmd[1] == "build":  # Record the genomes as SKA2 would read them (paths relative to its folder)
            for line in (kwargs["cwd"] / cmd[cmd.index("-f") + 1]).read_text().splitlines():
                name, path = line.split()
                seen[name] = len(next(read_records(kwargs["cwd"] / path)).seq)
        elif cmd[1] == "align":
            Path(cmd[cmd.index("-o") + 1]).write_text("".join(f">{n}\nA\n" for n in seen))

    monkeypatch.setattr(compare, "run", fake_run)
    ref = tmp_path / "ref.fasta"
    ref.write_text(">chr\n" + "ACGT" * 25 + "\n")
    asm = {}
    for name, header in [("a", "a_1 circular=true"), ("b", "b_1 circular=true"), ("pub", "pub")]:
        asm[name] = tmp_path / f"{name}.fasta"
        asm[name].write_text(f">{header}\n" + "ACGT" * 25 + "\n")
    (tmp_path / "logs").mkdir()
    compare.run_ska(ref, asm, tmp_path / "out", tmp_path / "logs", threads=1, min_freq=1.0, kmer=31,
                    added={"pub"})
    assert seen == {"Reference": 130, "a": 130, "b": 130, "pub": 130}  # All extended by 30 bases


def test_snp_alignment_keeps_only_columns_with_two_nucleotides(tmp_path):
    from bacon.compare import clean_alignment
    src = tmp_path / "raw.fasta"
    #                    SNP  ambiguity only  gap only  constant  SNP with N
    src.write_text(">Reference\nA" "A" "A" "A" "A\n>a\nC" "R" "-" "A" "N\n>b\nA" "A" "A" "A" "G\n")
    records = clean_alignment(src, tmp_path / "out.fasta", snps_only=True)
    assert [r.seq for r in records] == ["AA", "CN", "AG"]
    assert [r.seq for r in clean_alignment(src, tmp_path / "all.fasta")] == ["AAAAA", "CN-AN", "AAAAG"]


def test_ska_pan_alignment_has_only_snp_sites(tmp_path):
    """With the real SKA2 and --ska-min-freq below 1: a genome lacking a region does not make SNP sites."""
    import random
    import shutil

    import pytest
    if not shutil.which("ska"):
        pytest.skip("ska not installed")
    from bacon.compare import run_ska
    from bacon.seqio import read_records
    rng = random.Random(2)
    ref = "".join(rng.choice("ACGT") for _ in range(4000))
    (tmp_path / "ref.fasta").write_text(f">chr\n{ref}\n")
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    assemblies = {}
    for name in ("a", "b", "c", "d"):
        g = list(ref)
        if name == "a":
            g[1000] = swap[g[1000]]  # The only SNP
        if name == "b":
            g[3000] = "R"  # An ambiguous base: no SNP
        g = "".join(g)
        if name == "d":
            g = g[:2000] + g[2600:]  # Lacks 600 bp
        assemblies[name] = tmp_path / f"{name}.fasta"
        assemblies[name].write_text(f">{name}\n{g}\n")
    (tmp_path / "logs").mkdir()
    _, snps = run_ska(tmp_path / "ref.fasta", assemblies, tmp_path / "out", tmp_path / "logs", threads=1,
                      min_freq=0.5)
    assert {len(r.seq) for r in read_records(snps)} == {1}


@pytest.mark.parametrize("method", ["ska", "parsnp"])
def test_real_tools_with_spaces_and_accents_in_the_paths(tmp_path, method):
    """SKA2 splits its file list on whitespace; Parsnp its command lines (with the real programs)."""
    import random
    import shutil

    from bacon.compare import run_parsnp, run_ska, write_vcf
    if not all(shutil.which(p) for p in ({"ska": ["ska"], "parsnp": ["parsnp", "harvesttools"]}[method])):
        pytest.skip(f"{method} not installed")
    folder = tmp_path / "my r\u00e9sults"
    folder.mkdir()
    rng = random.Random(3)
    ref = "".join(rng.choice("ACGT") for _ in range(6000))
    (folder / "ref erence.fasta").write_text(f">chr\n{ref}\n")
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    assemblies = {}
    for i, name in enumerate(("a", "b", "c")):
        g = list(ref)
        g[1000 + 1000 * i] = swap[g[1000 + 1000 * i]]
        assemblies[name] = folder / f"gen ome {name}.fasta"
        assemblies[name].write_text(f">{name}\n{''.join(g)}\n")
    (folder / "lo gs").mkdir()
    out = folder / "com pared"
    if method == "ska":
        _, snps = run_ska(folder / "ref erence.fasta", assemblies, out, folder / "lo gs", threads=1, min_freq=1)
    else:
        _, snps = run_parsnp(folder / "ref erence.fasta", assemblies, out, folder / "lo gs", threads=1)
    assert {r.name: len(r.seq) for r in read_records(snps)} == {"Reference": 3, "a": 3, "b": 3, "c": 3}
    _, count = write_vcf(method, folder / "ref erence.fasta", out, folder / "lo gs", threads=1,
                         assemblies=assemblies, source="test")
    assert count == 3


def test_clean_vcf_rejects_positions_shifted_by_old_parsnp(tmp_path):
    """Parsnp 2.1.1 writes every SNP one base off: the reference then never carries its own allele."""
    shifted = "".join(f"chr\t{p}\t.\tA\tC,G\t40\tPASS\tNA\tGT\t1\t{g}\t1\n" for p, g in ((11, 2), (21, 0), (31, 1)))
    from bacon import BaconError
    with pytest.raises(BaconError, match="Parsnp 2.1.2 or later"):
        _clean(tmp_path, _vcf(tmp_path, shifted))
    one_ambiguous = "chr\t11\t.\tA\tC\t40\tPASS\tNA\tGT\t1\t0\t1\nchr\t21\t.\tA\tC\t40\tPASS\tNA\tGT\t0\t0\t1\n"
    assert _clean(tmp_path, _vcf(tmp_path, one_ambiguous))[0] == 1
