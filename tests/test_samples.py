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
    assert [f.name for f in samples["A"].files] == ["a1.fastq.gz", "a2.fastq.gz"]  # a2 listed twice: used once
    assert samples["B"].files == [tmp_path / "b.fastq.gz"]


def test_sample_sheet_paths_keep_their_inner_spaces(tmp_path, fastq):
    (tmp_path / "reads").mkdir()
    fastq("reads/a  b.fastq.gz", [("a", "ACGT")])
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tfile\nab \t reads/a  b.fastq.gz \n")
    [sample] = read_sample_sheet(sheet)
    assert sample.name == "ab" and sample.files == [tmp_path / "reads" / "a  b.fastq.gz"]


def test_sample_sheet_errors(tmp_path):
    sheet = tmp_path / "s.tsv"
    sheet.write_text("name\tpath\nA\tx.fq\n")
    with pytest.raises(BaconError, match="needs the columns"):
        read_sample_sheet(sheet)
    sheet.write_text("sample\tfile\nA\tmissing.fq\n")
    with pytest.raises(BaconError, match="file not found"):
        read_sample_sheet(sheet)


def test_files_in_hidden_subfolders_are_ignored(tmp_path, fastq):
    from bacon.samples import discover
    (tmp_path / "in" / "barcode01" / ".hidden").mkdir(parents=True)
    for path in ("in/barcode01/a.fastq", "in/barcode01/.hidden/b.fastq"):
        (tmp_path / path).write_text("@r\nACGT\n+\nIIII\n")
    [sample] = discover(tmp_path / "in")
    assert [f.name for f in sample.files] == ["a.fastq"]


def test_sheets_are_read_as_0_3_5_read_them(tmp_path):
    from bacon.metadata import sheet_metadata
    for name in ("a.fastq", "a b.fastq"):
        (tmp_path / name).write_text("@r\nACGT\n+\nIIII\n")
    sheet = tmp_path / "s.tsv"
    # Cells entirely in quotes (a spreadsheet's TSV export) lose them, in the header too, and "" inside is a quote;
    # a line starting with # after the header is a comment (a sample left out), for the sheet's metadata too
    sheet.write_text('"sample"\t"file"\t"note"\n"a"\t"a.fastq"\t""\n"b"\t"a b.fastq"\t"5"" tube, ""fine"""\n'
                     "#c\ta.fastq\tx\n")
    samples = read_sample_sheet(sheet)
    assert [(s.name, [f.name for f in s.files]) for s in samples] == [("a", ["a.fastq"]), ("b", ["a b.fastq"])]
    assert sheet_metadata(sheet).rows == {"a": {"note": ""}, "b": {"note": '5" tube, "fine"'}}
    sheet.write_text("sample\tfile\na\ta.fastq\n# b\ta b.fastq\n")
    assert [s.name for s in read_sample_sheet(sheet)] == ["a"]
    sheet.write_text("sample\tfile\tnote\tnote\na\ta.fastq\tx\ty\n")  # Still an error, saying what to do
    with pytest.raises(BaconError, match="column name 'note' is used more than once .* give each column"):
        read_sample_sheet(sheet)
