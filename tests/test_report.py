import json

from bacon.multiqc import distance_heatmap, reads_bargraph, sample_table
from bacon.report import _fasta_input, _tool, heatmap, identical_groups, methods_text


def test_identical_groups_are_transitive_and_sorted():
    names = ["a", "b", "c", "d", "e"]
    m = {x: {y: 0 if x == y else 5 for y in names} for x in names}
    for x, y in [("a", "b"), ("b", "c"), ("d", "e")]:
        m[x][y] = m[y][x] = 0
    assert identical_groups(names, m) == [["a", "b", "c"], ["d", "e"]]


def test_heatmap_cells_and_escaping():
    names = ["<a>", "b"]
    out = heatmap(names, {"<a>": {"<a>": 0, "b": 3}, "b": {"<a>": 3, "b": 0}})
    assert "&lt;a&gt;" in out and "<a>" not in out
    assert out.count("<td") == 4 and ">3</td>" in out and "3 SNPs" in out


def test_tool_version_number_only():
    info = {"tools": {"filtlong": {"version": "Filtlong v0.3.1"}, "minimap2": {"version": "2.31-r1302"},
                      "x": {"version": ""}}}
    assert _tool(info, "filtlong") == " 0.3.1"
    assert _tool(info, "minimap2") == " 2.31-r1302"
    assert _tool(info, "x") == "" and _tool(info, "missing") == ""


def test_methods_text_follows_the_settings():
    base = {"bacon_version": "0.3.0", "reference": {"file": "/x/ref.fasta", "length": 1000},
            "samples": {"s": {"files": ["s.fastq.gz"]}}, "comparison": {"tree": "t.nwk"}}
    text = methods_text({**base, "settings": {"baiting": "minimap2", "min_read_length": 500, "keep_percent": 95.0,
                                              "target_depth": 100, "assembler": "samtools", "snp_method": "ska",
                                              "ska_min_freq": 1.0, "tree": "fasttree"}})
    assert "samtools consensus" in text and "core SNPs" in text and "FastTree" in text
    assert "circular contigs" not in text and "fasta" not in text.lower().replace("ref.fasta", "")
    text = methods_text({**base, "samples": {"s": {"files": ["s.fa.gz"]}},
                         "settings": {"baiting": "bbduk", "kmer": 31, "min_read_length": 500, "keep_percent": 95.0,
                                      "target_depth": 100, "assembler": "flye", "read_type": "nano-hq",
                                      "flye_iterations": 3, "snp_method": "ska", "ska_min_freq": 0.5,
                                      "tree": "iqtree"}})
    assert "BBDuk" in text and "Flye" in text and "at least 50%" in text and "circular contigs" in text
    assert "IQ-TREE" in text and "without qualities" in text


def test_fasta_input_detection():
    assert not _fasta_input({"samples": {"a": {"files": ["a.fastq.gz", "b.fq"]}}})
    assert _fasta_input({"samples": {"a": {"files": ["a.fna"]}}})


ROWS = [
    {"Sample": "s1", "Status": "ok", "Raw_bases": "1000", "Baited_bases": "600", "Filtered_bases": "500",
     "Baited_reads": "6", "Est_depth": "12.5", "Note": "low depth (12x)", "Contigs": "1"},
    {"Sample": "s2", "Status": "failed (bait)", "Raw_bases": "NA", "Baited_bases": "NA", "Filtered_bases": "NA",
     "Note": "no reads matched the reference"},
]


def test_multiqc_table_and_bargraph():
    table = sample_table(ROWS)
    assert table["plot_type"] == "table" and table["data"]["s1"]["Est_depth"] == 12.5
    assert table["data"]["s1"]["Contigs"] == 1 and "Raw_reads" not in table["data"]["s1"]
    assert table["data"]["s2"] == {"Status": "failed (bait)", "Note": "no reads matched the reference"}
    bars = reads_bargraph(ROWS)
    assert bars["data"] == {"s1": {"Kept": 500, "Baited, filtered out": 100, "Off-target": 400}}
    assert reads_bargraph(ROWS[1:]) is None
    json.dumps(table)


def test_multiqc_heatmap_in_tree_order(tmp_path):
    dist = tmp_path / "d.tsv"
    dist.write_text("snp-dists\tReference\ta\tb\nReference\t0\t1\t5\na\t1\t0\t4\nb\t5\t4\t0\n")
    tree = tmp_path / "t.nwk"
    tree.write_text("(b:1,(a:1,Reference:1):1);\n")
    hm = distance_heatmap(dist, tree)
    assert hm["xcats"] == ["b", "a", "Reference"] == hm["ycats"]
    assert hm["data"] == [[0, 4, 5], [4, 0, 1], [5, 1, 0]]
    assert distance_heatmap(dist)["xcats"] == ["Reference", "a", "b"]
