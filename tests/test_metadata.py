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
    sort_key,
    unusable_reason,
    write_copy,
)


def test_read_table_formats(tmp_path):
    tsv = tmp_path / "m.tsv"
    tsv.write_bytes(b"\xef\xbb\xbfSample\tGroup\t\r\n# a comment\r\n\r\n s1 \t A \r\ns2\tB\textra\r\n")
    header, rows = read_table(tsv, "Metadata file")
    assert header == ["Sample", "Group"]  # BOM gone, names stripped, the empty trailing column dropped
    assert rows == [{"Sample": "s1", "Group": "A", "": ""}, {"Sample": "s2", "Group": "B", "": "extra"}]
    csv = tmp_path / "m.csv"
    csv.write_text('sample,note\ns1,"a, quoted\tvalue"\n')
    assert read_table(csv, "x")[1] == [{"sample": "s1", "note": "a, quoted value"}]  # Tabs become spaces
    with pytest.raises(BaconError, match="Metadata file not found"):
        read_table(tmp_path / "none.tsv", "Metadata file")
    (tmp_path / "empty.tsv").write_text("# only a comment\n\n")
    with pytest.raises(BaconError, match="is empty"):
        read_table(tmp_path / "empty.tsv", "Metadata file")


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
    sheet = Metadata(["site", "group"], {"s1": {"site": "n", "group": "X"}, "s3": {"site": "s", "group": "Y"}},
                     ["w"])
    m = merge(given, sheet)
    assert m.columns == ["group", "year", "site"]
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
    assert is_numeric(["1", "", "2.5"]) and not is_numeric(["1", "x"]) and not is_numeric(["", ""])
    assert not is_numeric(["nan"])


def test_write_copy(tmp_path):
    m = Metadata(["g", "n"], {"a": {"g": "A", "n": ""}, "b": {"g": "", "n": "2"}})
    path = tmp_path / "metadata.tsv"
    write_copy(path, m)
    assert path.read_text() == "sample\tg\tn\na\tA\t\nb\t\t2\n"
    before = path.stat().st_mtime_ns
    write_copy(path, m)
    assert path.stat().st_mtime_ns == before  # Unchanged: not rewritten
    assert read_metadata(path).rows == m.rows  # The copy reads back as a metadata file
