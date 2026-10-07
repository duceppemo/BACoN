"""Sample metadata: parsing, sample-sheet columns, merging, the colour column."""

import pytest

from bacon import BaconError
from bacon.metadata import (
    Metadata,
    choose_colour_column,
    is_numeric,
    merge,
    read_metadata,
    read_table,
    restrict,
    sheet_metadata,
    shown_name,
    sort_key,
    unusable_reason,
    write_copy,
)


def test_read_table_formats(tmp_path):
    tsv = tmp_path / "m.tsv"
    tsv.write_bytes(b"\xef\xbb\xbf# a comment\r\n\r\nSample\tGroup\t\r\n\r\n s1 \t A \r\ns2\tB\textra\r\n")
    header, rows = read_table(tsv, "Metadata file")
    assert header == ["Sample", "Group"]  # BOM gone, names stripped, the empty trailing column dropped
    assert rows == [{"Sample": "s1", "Group": "A"}, {"Sample": "s2", "Group": "B"}]  # Extra cells ignored
    csv = tmp_path / "m.csv"
    csv.write_text('sample,note\ns1,"a, quoted\tvalue"\n')
    assert read_table(csv, "x")[1] == [{"sample": "s1", "note": "a, quoted\tvalue"}]  # Quotes: commas kept
    with pytest.raises(BaconError, match="Metadata file not found"):
        read_table(tmp_path / "none.tsv", "Metadata file")
    (tmp_path / "empty.tsv").write_text("# only a comment\n\n")
    with pytest.raises(BaconError, match="is empty"):
        read_table(tmp_path / "empty.tsv", "Metadata file")


def test_read_table_refuses_unclosed_quotes_and_quoted_line_breaks(tmp_path):
    csv = tmp_path / "m.csv"
    csv.write_text('sample,size,group\ns1,"5 inch,A\ns2,3,B\ns3,2,B\n')  # Would swallow s2 and s3 silently
    with pytest.raises(BaconError, match=r"m.csv, line 2: .*unclosed quote"):
        read_table(csv, "Metadata file")
    csv.write_text('sample,note\ns1,"line one\nline two"\ns2,x\n')  # A line break in a value breaks the copy
    with pytest.raises(BaconError, match="line 2: a quoted value spans several lines"):
        read_table(csv, "Metadata file")
    csv.write_text('\n# c\nsample,note\ns1,"say ""hi"""\n\ns2,"ab"c\n')  # The line is the file's
    with pytest.raises(BaconError, match="line 6: "):
        read_table(csv, "Metadata file")
    tsv = tmp_path / "m.tsv"
    tsv.write_text('sample\tsize\tgroup\ns1\t"5 inch\tA\ns2\t3\tB\ns3\t2\tB\n')
    assert [r["sample"] for r in read_table(tsv, "x")[1]] == ["s1", "s2", "s3"]  # In a TSV a tab always splits
    assert read_table(tsv, "x")[1][0]["size"] == '"5 inch'  # Not entirely quoted: kept as written


def test_tsv_cells_entirely_in_quotes_lose_them(tmp_path):
    # As a spreadsheet exports them, and as 0.3.5 read them (through the csv module); "" inside is one quote
    tsv = tmp_path / "m.tsv"
    tsv.write_text('"sample"\t"colour"\t"note"\n"s1"\t"#FF0000"\t"a ""b"" c"\ns2\t""\t5" tube\n')
    header, rows = read_table(tsv, "x")
    assert header == ["sample", "colour", "note"]
    assert rows == [{"sample": "s1", "colour": "#FF0000", "note": 'a "b" c'},
                    {"sample": "s2", "colour": "", "note": '5" tube'}]


def test_read_table_strips_and_only_metadata_values_lose_inner_whitespace(tmp_path):
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample\tfile\tsite\nab \t reads/a  b.fastq.gz\tnorth   shore\n")
    header, rows = read_table(sheet, "Sample sheet")
    assert rows == [{"sample": "ab", "file": "reads/a  b.fastq.gz", "site": "north   shore"}]  # Paths kept whole
    assert sheet_metadata(sheet).rows == {"ab": {"site": "north shore"}}  # Metadata values: one space
    csv = tmp_path / "m.csv"
    csv.write_text("sample,no   te\ns1,a\tb\n")
    assert read_metadata(csv).rows == {"s1": {"no te": "a b"}}  # A tab in a CSV cell would break the copy


def test_comment_lines_count_only_before_the_header(tmp_path):
    path = tmp_path / "m.csv"
    path.write_text("# simulated values\n\n# second comment\ncolour,sample\n#FF0000,s1\n#00FF00,s2\n")
    assert read_metadata(path).rows == {"s1": {"colour": "#FF0000"}, "s2": {"colour": "#00FF00"}}
    path.write_text("sample,g\ns1,A\n# not a comment\n")
    assert list(read_metadata(path).rows) == ["s1", "# not a comment"]
    assert [r["sample"] for r in read_table(path, "Sample sheet", comment_lines=True)[1]] == ["s1"]  # A sheet


def test_duplicate_column_names_are_an_error(tmp_path):
    path = tmp_path / "m.csv"
    path.write_text("sample,Group,group,g\ns1,A,B,C\n")
    with pytest.raises(BaconError, match="column name 'Group', 'group' is used more than once .*its own name"):
        read_metadata(path)
    path.write_text("sample,g,Sample\ns1,A,B\n")  # A second 'sample' column, whatever its case
    with pytest.raises(BaconError, match="column name 'sample', 'Sample' is used more than once"):
        read_metadata(path)
    path.write_text("sample,file,File\ns1,a.fq,b.fq\n")
    with pytest.raises(BaconError, match="Sample sheet .*column name 'file', 'File' is used more than once"):
        sheet_metadata(path)
    path.write_text("sample,g,,\ns1,A,,\n")  # Several empty names are not duplicates
    assert read_metadata(path).columns == ["g"]


def test_shown_name_marks_the_tables_own_columns():
    assert shown_name("Status", ["Sample", "Status", "Note"]) == "Status (metadata)"
    assert shown_name("note", ["Sample", "Status", "Note"]) == "note (metadata)"
    assert shown_name("Site", ["Sample", "Status", "Note"]) == "Site"


def test_read_metadata_values_and_warnings(tmp_path):
    path = tmp_path / "m.tsv"
    path.write_text("SAMPLE\tgroup\tyear\ts1\tA\t2020\ns2\tNA\t-\ns3\tna\t\ns1\tB\t1999\n\tC\t1\n".replace(
        "year\ts1", "year\ns1"))
    m = read_metadata(path)
    assert m.columns == ["group", "year"]
    assert m.rows == {"s1": {"group": "A", "year": "2020"}, "s2": {"group": "", "year": ""},
                      "s3": {"group": "", "year": ""}}  # First row of s1 kept; NA, na, - and "" are missing
    assert len(m.warnings) == 1 and "listed more than once" in m.warnings[0] and "s1" in m.warnings[0]
    assert m.value("s1", "group") == "A" and m.value("nobody", "group") == "" and m.values("year") == ["2020", "", ""]


def test_read_metadata_errors(tmp_path):
    path = tmp_path / "m.csv"
    path.write_text("name,group\ns1,A\n")
    with pytest.raises(BaconError, match="needs a 'sample' column"):
        read_metadata(path)
    path.write_text("sample\ns1\n")
    with pytest.raises(BaconError, match="no column besides"):
        read_metadata(path)


def test_sheet_metadata_agreeing_and_conflicting_rows(tmp_path):
    sheet = tmp_path / "s.csv"
    sheet.write_text("Sample,File,site,year\nA,a1.fq,north,2020\nA,a2.fq,north,\nB,b.fq,south,2021\n")
    m = sheet_metadata(sheet)
    assert m.columns == ["site", "year"] and m.rows == {"A": {"site": "north", "year": "2020"},
                                                        "B": {"site": "south", "year": "2021"}}
    assert m.warnings == []
    sheet.write_text("sample\tfile\tsite\nA\ta1.fq\tnorth\nA\ta2.fq\tsouth\n")
    m = sheet_metadata(sheet)
    assert m.rows["A"]["site"] == "north" and "different values" in m.warnings[0] and "'north' kept" in m.warnings[0]
    sheet.write_text("sample\tfile\nA\ta1.fq\n")
    assert sheet_metadata(sheet) is None


def test_merge_gives_the_first_table_precedence_column_by_column():
    given = Metadata(["group", "year"], {"s1": {"group": "A", "year": ""}, "s2": {"group": "B", "year": "2"}})
    sheet = Metadata(["site", "Group"], {"s1": {"site": "n", "Group": "X"}, "s3": {"site": "s", "Group": "Y"}},
                     ["w"])  # 'Group' is the sheet's name for the same column
    m = merge(given, sheet)
    assert m.columns == ["group", "year", "site"]  # Not 'Group' again
    assert m.rows == {"s1": {"group": "A", "year": "", "site": "n"}, "s2": {"group": "B", "year": "2", "site": ""},
                      "s3": {"group": "", "year": "", "site": "s"}}  # s3's group comes from --metadata: none
    assert m.warnings == ["w"]
    assert merge(given, None) is given and merge(None, sheet) is sheet and merge(None, None) is None


def test_restrict_to_the_samples_of_the_run():
    m = Metadata(["g"], {"b": {"g": "1"}, "x": {"g": "2"}, "a": {"g": "3"}})
    r, unmatched = restrict(m, ["a", "b", "c"])
    assert list(r.rows) == ["a", "b", "c"] and r.rows["c"] == {"g": ""} and unmatched == ["x"]


def test_unusable_reason():
    assert unusable_reason(["A", "B", "", "A"]) is None
    assert unusable_reason(["", "", ""]) == "no value"
    assert "9 distinct" in unusable_reason([str(i) for i in range(9)])
    assert unusable_reason(["A"] * 9 + ["B"]) is None
    assert unusable_reason([str(i) for i in range(8)]) is None  # 8 distinct values: the palette's size
    assert "free text" in unusable_reason([f"v{i % 6}" for i in range(10)])  # 6 distinct among 10
    assert unusable_reason([f"v{i % 5}" for i in range(10)]) is None  # 5 of 10 is still a category
    assert unusable_reason([f"v{i}" for i in range(6)]) is None  # Below 10 samples, up to 8 distinct is fine
    assert "longer than 30" in unusable_reason(["a sentence that goes on and on and on", "and another one that is as long"])
    assert unusable_reason(["a sentence that goes on and on and on", "short"]) is None  # 21 characters on average


def test_choose_colour_column():
    m = Metadata(["comment", "Group", "n"], {f"s{i}": {"comment": f"free text number {i} about this sample",
                                                       "Group": "A" if i % 2 else "B", "n": str(i)}
                                             for i in range(12)})
    assert choose_colour_column(m, None) == ("Group", None)
    two = Metadata(["Group", "site"], {k: {"Group": v["Group"], "site": "x" if k < "s5" else "y"}
                                       for k, v in m.rows.items()})
    assert choose_colour_column(two, None) == ("Group", None)  # The first usable column, not the last
    assert choose_colour_column(m, "group") == ("Group", None)  # Any case
    assert choose_colour_column(m, "NONE") == (None, None)
    column, why = choose_colour_column(m, "n")
    assert column is None and "'n' cannot colour" in why and "12 distinct" in why
    with pytest.raises(BaconError, match="no such metadata column"):
        choose_colour_column(m, "missing")
    only_text = Metadata(["comment"], {k: {"comment": v["comment"]} for k, v in m.rows.items()})
    column, why = choose_colour_column(only_text, None)
    assert column is None and why.startswith("no metadata column can colour")


def test_sort_key_and_is_numeric():
    assert sorted(["10", "9", "b", "A", "2.5"], key=sort_key) == ["2.5", "9", "10", "A", "b"]
    assert sorted(["v10", "v9", "V2", "a10b", "a9b"], key=sort_key) == ["a9b", "a10b", "V2", "v9", "v10"]  # Natural
    assert is_numeric(["1", "", "2.5"]) and not is_numeric(["1", "x"]) and not is_numeric(["", ""])
    assert is_numeric(["1e5", "0012", ".5", "+5", "-1.", "2E-3"])  # What the page's parseFloat reads whole
    for value in ("nan", "inf", "Infinity", "1_000", "1,000", "0x10", "١٢", "1e", "e5"):
        assert not is_numeric([value]), value
        assert sort_key(value)[0] == 1, value  # Sorted as text


def test_write_copy(tmp_path):
    m = Metadata(["g", "n"], {"a": {"g": "A", "n": ""}, "b": {"g": "", "n": "2"}})
    path = tmp_path / "metadata.tsv"
    write_copy(path, m)
    assert path.read_text() == "sample\tg\tn\na\tA\t\nb\t\t2\n"
    before = path.stat().st_mtime_ns
    write_copy(path, m)
    assert path.stat().st_mtime_ns == before  # Unchanged: not rewritten
    assert read_metadata(path).rows == m.rows  # The copy reads back as a metadata file
