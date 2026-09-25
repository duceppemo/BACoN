import pytest

from bacon import BaconError
from bacon.samples import discover, read_sample_sheet


def test_single_file(fastq):
    samples = discover(fastq("s1.fastq.gz", [("r", "ACGT")]))
    assert [(s.name, s.fmt) for s in samples] == [("s1", "fastq")]


def test_folder_files_and_subfolders(tmp_path, fastq, fasta):
    d = tmp_path / "in"
    (d / "barcode01").mkdir(parents=True)
    (d / "unclassified").mkdir()
    (d / "notes").mkdir()
    (d / "notes" / "readme.txt").write_text("hi")
    fastq("in/barcode01/chunk_0.fastq.gz", [("a", "ACGT")])
    fastq("in/barcode01/chunk_1.fastq.gz", [("b", "ACGT")])
    fastq("in/unclassified/x.fastq.gz", [("c", "ACGT")])
    fasta("in/asm.fasta", [("c1", "ACGT")])
    (d / "README.md").write_text("ignored")
    samples = {s.name: s for s in discover(d)}
    assert sorted(samples) == ["asm", "barcode01"]
    assert len(samples["barcode01"].files) == 2
    assert samples["asm"].fmt == "fasta"


def test_duplicate_names_rejected(tmp_path, fastq):
    (tmp_path / "in").mkdir()
    fastq("in/s1.fastq.gz", [("a", "ACGT")])
    fastq("in/s1.fq", [("a", "ACGT")], gz=False)
    with pytest.raises(BaconError, match="same sample name"):
        discover(tmp_path / "in")


def test_mixed_formats_rejected(tmp_path, fastq, fasta):
    (tmp_path / "in" / "s").mkdir(parents=True)
    fastq("in/s/a.fastq.gz", [("a", "ACGT")])
    fasta("in/s/b.fasta", [("b", "ACGT")])
    with pytest.raises(BaconError, match="mixes fasta and fastq"):
        discover(tmp_path / "in")


def test_empty_files_skipped_and_nothing_left_is_an_error(tmp_path):
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "empty.fastq").write_text("")
    with pytest.raises(BaconError, match="No fasta/fastq"):
        discover(tmp_path / "in")


@pytest.mark.parametrize("name", ["bad name", "Reference", "all_assemblies"])
def test_invalid_or_reserved_names(tmp_path, fastq, name):
    fastq(f"{name}.fastq.gz", [("a", "ACGT")])
    with pytest.raises(BaconError, match="Invalid sample name|reserved"):
        discover(tmp_path / f"{name}.fastq.gz")


def test_sample_sheet(tmp_path, fastq):
    fastq("a1.fastq.gz", [("a", "ACGT")])
    fastq("a2.fastq.gz", [("b", "ACGT")])
    fastq("b.fastq.gz", [("c", "ACGT")])
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("# comment\nSample,File\nA,a1.fastq.gz;a2.fastq.gz\nB,b.fastq.gz\nA,a2.fastq.gz\n")
    samples = {s.name: s for s in read_sample_sheet(sheet)}
    assert [f.name for f in samples["A"].files] == ["a1.fastq.gz", "a2.fastq.gz", "a2.fastq.gz"]
    assert samples["B"].files == [tmp_path / "b.fastq.gz"]


def test_sample_sheet_errors(tmp_path):
    sheet = tmp_path / "s.tsv"
    sheet.write_text("name\tpath\nA\tx.fq\n")
    with pytest.raises(BaconError, match="needs the columns"):
        read_sample_sheet(sheet)
    sheet.write_text("sample\tfile\nA\tmissing.fq\n")
    with pytest.raises(BaconError, match="file not found"):
        read_sample_sheet(sheet)
