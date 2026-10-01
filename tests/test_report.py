import json

from bacon.multiqc import distance_heatmap, reads_bargraph, sample_table
from bacon.report import (
    _fasta_samples,
    _tool,
    build_report,
    distinct_count,
    heatmap,
    identical_groups,
    methods_text,
)


def test_identical_groups_require_zero_distance_between_all_members():
    names = ["a", "b", "c", "d", "e"]
    m = {x: {y: 0 if x == y else 5 for y in names} for x in names}
    for x, y in [("a", "b"), ("b", "c"), ("a", "c"), ("d", "e")]:
        m[x][y] = m[y][x] = 0
    assert identical_groups(names, m) == [["a", "b", "c"], ["d", "e"]]
    assert distinct_count(names, m) == 2


def test_zero_distances_are_not_transitive():
    # b has N where a and c differ: a = b = c by distance, but a and c differ.
    m = {"a": {"a": 0, "b": 0, "c": 5}, "b": {"a": 0, "b": 0, "c": 0}, "c": {"a": 5, "b": 0, "c": 0}}
    groups = identical_groups(["a", "b", "c"], m)
    assert groups == [["a", "b"]] and not any("a" in g and "c" in g for g in groups)
    assert distinct_count(["a", "b", "c"], m) == 2


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


SETTINGS = {"baiting": "minimap2", "min_read_length": 500, "keep_percent": 95.0, "target_depth": 100,
            "assembler": "samtools", "snp_method": "ska", "ska_min_freq": 1.0, "tree": "fasttree",
            "template_gaps": "n", "read_type": "nano-hq", "flye_iterations": 3, "kmer": 31}
BASE = {"bacon_version": "0.3.1", "reference": {"file": "/x/ref.fasta", "length": 1000},
        "samples": {"s": {"files": ["s.fastq.gz"]}},
        "comparison": {"tree": "t.nwk", "distances": "d.tsv", "method": "ska", "core_snps": 3}}


def _methods(rows=None, **changes):
    info = {**BASE, "settings": {**SETTINGS, **changes.pop("settings", {})}, **changes}
    return methods_text(info, rows)


def test_methods_text_default():
    text = _methods()
    assert "samtools consensus 1" not in text  # No version recorded: no number invented
    assert "samtools consensus" in text and "core SNPs" in text and "FastTree" in text
    assert "best 95% were kept" in text and "Filtlong" in text
    assert "circular contigs" not in text and "fasta" not in text.lower().replace("ref.fasta", "")
    assert "SNP alignment" in text and "Pairwise SNP distances count" in text and "VCF" not in text
    assert "ska map" in _methods(comparison={**BASE["comparison"], "vcf": "snps.vcf"})


def test_methods_text_settings_combinations():
    text = _methods(samples={"s": {"files": ["s.fa.gz"]}}, settings={"baiting": "bbduk", "assembler": "flye",
                    "ska_min_freq": 0.5, "tree": "iqtree", "min_size": 3000, "genome_size": 150000})
    assert "BBDuk" in text and "Flye" in text and "at least 50% of the genomes" in text
    assert "best 95%" not in text and "filtered by BACoN" in text  # Fasta reads: no Filtlong, no quality cut
    assert "minimum overlap 3000 bp" in text and "genome size as given" in text
    assert "bootstrap consensus tree" in text and "circular contigs" not in text  # None circular
    text = _methods(rows=[{"Circular_contigs": "1"}], settings={"assembler": "myloasm", "ska_min_freq": 0.999})
    assert "circular contigs (and the reference" in text and "99.9%" in text
    text = _methods(samples={"a": {"files": ["a.fq"]}, "b": {"files": ["b.fasta"]}})
    assert "For samples with base qualities" in text and "without qualities" in text
    assert "replaced by the reference" in _methods(settings={"template_gaps": "reference"})
    text = _methods(settings={"snp_method": "parsnp", "add_genomes": ["x.fa", "y.fa"]})
    assert "Parsnp" in text and "core-genome alignment" in text and "2 finished genomes were added" in text


def test_methods_text_failed_or_skipped_comparison():
    for comparison in ({"failed": "ska build failed"}, {"failed": "FastTree failed", "distances": "d.tsv"},
                       {"skipped": "only 2 assemblies"}, None):
        text = _methods(comparison=comparison)
        assert "SKA2" not in text and "A tree was built" not in text


def test_fasta_samples():
    assert _fasta_samples({"samples": {"a": {"files": ["a.fastq.gz", "b.fq"]}}}) == (0, 1)
    assert _fasta_samples({"samples": {"a": {"files": ["a.fna"]}, "b": {"files": ["b.fq"]}}}) == (1, 2)


def _folder(tmp_path, comparison, summary="Sample\tStatus\na\tok\nb\tok\n", matrix=None, **info):
    out = tmp_path / "out"
    (out / "4_compared" / "ska").mkdir(parents=True)
    (out / "summary.tsv").write_text(summary)
    matrix = matrix or "snp-dists\tReference\ta\tb\nReference\t0\t2\t3\na\t2\t0\t1\nb\t3\t1\t0\n"
    (out / "4_compared" / "ska" / "snp_distances.tsv").write_text(matrix)
    (out / "run_info.json").write_text(json.dumps({"settings": SETTINGS, "comparison": comparison, **info}))
    return out


def test_report_zero_sites_claims_no_identity(tmp_path):
    zeros = "snp-dists\tReference\ta\tb\nReference\t0\t0\t0\na\t0\t0\t0\nb\t0\t0\t0\n"
    out = _folder(tmp_path, {"method": "ska", "core_snps": 0, "tree": None,
                             "distances": "4_compared/ska/snp_distances.tsv"}, matrix=zeros)
    page = build_report(out)
    assert "No SNP site was compared" in page and "Identical genomes" not in page and "distinct" not in page
    assert "see --ska-min-freq" in page and "midpoint-rooted" not in page


def test_report_finds_files_of_a_moved_folder(tmp_path):
    out = _folder(tmp_path, {"method": "ska", "core_snps": 3,
                             "distances": "/old/place/out/4_compared/ska/snp_distances.tsv"})
    page = build_report(out)
    assert "SNP distances" in page and "not found" not in page
    (out / "4_compared" / "ska" / "snp_distances.tsv").unlink()
    assert "/old/place/out/4_compared/ska/snp_distances.tsv was not found" in build_report(out)


def test_report_survives_incomplete_run_info(tmp_path):
    out = _folder(tmp_path, {"method": "ska", "distances": "4_compared/ska/snp_distances.tsv"},
                  summary="", tools={"x": None}, command_line=None)
    info = json.loads((out / "run_info.json").read_text())
    info["settings"] = {"keep_percent": None}
    (out / "run_info.json").write_text(json.dumps(info))
    page = build_report(out)
    assert "BACoN" in page and "in the order of the distance table" in page


ROWS = [
    {"Sample": "s1", "Status": "ok", "Raw_bases": "1000", "Baited_bases": "600", "Filtered_bases": "500",
     "Baited_reads": "6", "Est_depth": "12.5", "Note": "low depth (12x)", "Contigs": "1"},
    {"Sample": "s2", "Status": "failed (bait)", "Raw_bases": "NA", "Baited_bases": "NA", "Filtered_bases": "NA",
     "Note": "no reads matched the reference"},
]


def test_multiqc_table_and_bargraph():
    table = sample_table(ROWS, "run1")
    assert table["plot_type"] == "table" and table["data"]["s1"]["Est_depth"] == 12.5
    assert table["id"] == "bacon_samples_run1" and table["section_name"] == "BACoN run1: samples"
    assert table["data"]["s1"]["Contigs"] == 1 and "Raw_reads" not in table["data"]["s1"]
    assert table["data"]["s2"] == {"Status": "failed (bait)", "Note": "no reads matched the reference"}
    bars = reads_bargraph(ROWS, "run1")
    assert bars["data"] == {"s1": {"Kept": 500, "Baited, filtered out": 100, "Off-target": 400}}
    assert bars["id"] == "bacon_reads_run1"
    assert reads_bargraph(ROWS[1:]) is None
    bbduk = [{"Sample": "s3", "Raw_bases": "NA", "Baited_bases": "700", "Filtered_bases": "NA"}]
    assert reads_bargraph(bbduk)["data"] == {"s3": {"Baited": 700}}
    json.dumps(table)


def test_multiqc_heatmap_in_tree_order(tmp_path):
    dist = tmp_path / "d.tsv"
    dist.write_text("snp-dists\tReference\ta\tb\nReference\t0\t1\t5\na\t1\t0\t4\nb\t5\t4\t0\n")
    tree = tmp_path / "t.nwk"
    tree.write_text("(b:1,(a:1,Reference:1):1);\n")
    hm = distance_heatmap(dist, tree, "my run")
    assert hm["xcats"] == ["b", "a", "Reference"] == hm["ycats"]
    assert hm["data"] == [[0, 4, 5], [4, 0, 1], [5, 1, 0]]
    assert hm["id"] == "bacon_distances_my_run"
    assert distance_heatmap(dist)["xcats"] == ["Reference", "a", "b"]
    tree.write_text("(b:1,(a_renamed:1,Reference:1):1);\n")  # Names that do not match: table order
    assert distance_heatmap(dist, tree)["xcats"] == ["Reference", "a", "b"]


def test_multiqc_ids_differ_between_runs_and_stale_files_go(tmp_path):
    from bacon.multiqc import write_multiqc
    dist = tmp_path / "d.tsv"
    dist.write_text("snp-dists\tReference\ta\nReference\t0\t1\na\t1\t0\n")
    ids = set()
    for run in ("run A", "run-B"):
        out = tmp_path / run
        out.mkdir()
        write_multiqc(out, ROWS, dist)
        ids |= {json.loads(p.read_text())["id"] for p in out.glob("*_mqc.json")}
    assert len(ids) == 6
    write_multiqc(tmp_path / "run A", ROWS, None)
    assert not (tmp_path / "run A" / "bacon_distances_mqc.json").exists()


def test_no_inline_heatmap_for_hundreds_of_genomes(tmp_path):
    names = ["Reference"] + [f"g{i}" for i in range(160)]
    rows = ["\t".join([n] + ["1" if n != m else "0" for m in names]) for n in names]
    matrix = "\t".join(["snp-dists", *names]) + "\n" + "\n".join(rows) + "\n"
    out = _folder(tmp_path, {"method": "ska", "core_snps": 5, "distances": "4_compared/ska/snp_distances.tsv"},
                  matrix=matrix)
    page = build_report(out)
    assert 'class="heat"' not in page and "too many for a heatmap" in page and "4_compared/ska/snp_distances.tsv" in page


def test_methods_text_gives_bbduk_mismatches():
    from bacon.report import methods_text
    for hdist, text in [(0, "no mismatch"), (1, "up to one mismatch"), (2, "up to 2 mismatches"),
                        (None, "up to 2 mismatches")]:
        settings = {"baiting": "bbduk", "kmer": 31, **({"hdist": hdist} if hdist is not None else {})}
        assert f"31-mer ({text})" in methods_text({"settings": settings})
