import json
import re

from bacon.multiqc import distance_heatmap, reads_bargraph, sample_table
from bacon.newick import parse
from bacon.report import (
    MAX_LABELLED_CELLS,
    Bar,
    _fasta_samples,
    _tool,
    bar_chart,
    build_report,
    depth_bars,
    distance_bins,
    distinct_count,
    genome_map,
    heatmap,
    identical_groups,
    methods_text,
    n_bars,
    n_per_bin,
    read_vcf,
    tree_svg,
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
    assert out.count("<title>&lt;a&gt; – b: 3 SNPs</title>") == 1 and out.count("</title></rect>") == 4
    assert ">3</text>" in out and "1 SNP</title>" not in out


def test_distance_bins_adapt_to_the_range():
    assert distance_bins(0) == [(0, 0)]
    assert distance_bins(1) == [(0, 0), (1, 1)]
    assert distance_bins(2) == [(0, 0), (1, 1), (2, 2)]
    bins = distance_bins(118)
    assert bins[0] == (0, 0) and bins[-1][1] == 118 and len(bins) == 5
    assert all(b[0] == a[1] + 1 for a, b in zip(bins, bins[1:]))  # Contiguous
    assert [hi for _, hi in bins] == sorted(hi for _, hi in bins)


def test_heatmap_groups_bins_and_legend():
    names = ["a", "b", "c", "d"]
    m = {x: {y: 0 if x == y else 12 for y in names} for x in names}
    m["a"]["b"] = m["b"]["a"] = 0
    groups = identical_groups(names, m)
    out = heatmap(names, m, groups)
    assert "group 1 (2)" in out and 'fill="var(--s1)"' in out  # The block on the right and the bands
    assert out.count("<title>a: group 1</title>") == 2 and "c: group" not in out  # Bands on both axes; c alone
    assert "SNPs:" in out and ">0</text>" in out and ">12</text>" in out  # Legend labels
    assert 'class="q0"' in out and 'class="q4"' in out  # Zero and the largest class
    assert "<title>a – c: 12 SNPs</title>" in out
    without = heatmap(names, m)
    assert "group" not in without and "var(--s" not in without


def test_heatmap_values_on_hover_only_for_many_genomes():
    names = [f"g{i}" for i in range(MAX_LABELLED_CELLS + 1)]
    m = {x: {y: 0 if x == y else 1 for y in names} for x in names}
    out = heatmap(names, m)
    assert "<title>g0 – g1: 1 SNP</title>" in out
    assert out.count("</text>") == 2 * len(names) + 1 + 2  # Row and column labels, "SNPs:" and two legend classes
    assert '">1</text>' in heatmap(names[:3], {x: m[x] for x in names[:3]})


def test_tree_svg_marks_groups_supports_and_scale():
    root = parse("((a:0.01,b:0.01)0.95:0.02,(Reference:0.005,c:0.005)0.80:0.03);")
    out = tree_svg(root, {"a": 0, "b": 0}, "ref.fa", 100)
    assert out.count('fill="var(--s1)"') == 2 and "<title>group 1</title>" in out
    assert 'font-weight="600">Reference <tspan' in out and "ref.fa" in out
    assert ">0.95</text>" in out and ">0.80</text>" in out
    assert "substitutions per site (about 1 SNP)" in out  # 0.01 per site x 100 sites
    assert "(about" not in tree_svg(root, {}, "", None)
    flat = tree_svg(parse("(a:0,b:0);"), {}, "", 10)
    assert "substitutions" not in flat  # No scale bar when every branch is 0


def test_tree_svg_escapes_names():
    out = tree_svg(parse("('<x>':0.1,Reference:0.2)'<s>':0.1;"), {"<x>": 0}, "<r>", 10)
    assert "&lt;x&gt;" in out and "&lt;r&gt;" in out and "<x>" not in out and "<r>" not in out
    assert "(about 0.5 SNPs)" in out  # 0.05 per site x 10 sites


def test_bar_chart_sorts_flags_and_lists_failed_samples():
    rows = [
        {"Sample": "low", "Status": "ok", "Est_depth": "12.5", "Filtered_reads": "10", "Filtered_N50": "500",
         "Note": "low depth (12x)", "N_bases": "7", "Assembly_length": "1000"},
        {"Sample": "high", "Status": "ok", "Est_depth": "80", "N_bases": "0", "Assembly_length": "1000"},
        {"Sample": "gone", "Status": "failed (bait)", "Est_depth": "NA", "N_bases": "NA",
         "Note": "no reads matched"},
        {"Sample": "<odd>", "Status": "ok", "Est_depth": "NA", "N_bases": "NA"},
    ]
    out = bar_chart(depth_bars(rows), "Depth", unit="x", threshold=20, threshold_text="20x flag")
    assert out.index(">high<") < out.index(">low<") < out.index(">&lt;odd&gt;<") < out.index(">gone<")
    assert out.count('class="bar"') == 2 and 't-ink">12.5x</text>' in out and 't-ink">80x</text>' not in out
    assert 'class="thresh"' in out and "20x flag" in out
    assert "failed (bait)" in out and "no reads matched" in out and ">NA<" in out
    assert "low depth (12x)" in out and "10 reads, N50 500" in out and "<odd>" not in out
    out = bar_chart(n_bars(rows), "N bases", integers=True)
    assert ">7</text>" in out and "7 N bases of 1,000 bp (0.70%)" in out and 'class="bar0"' in out
    assert ">0.2<" not in out  # Integer ticks
    assert bar_chart([], "Empty").count("<svg") == 1


def test_bar_chart_empty_and_threshold_only():
    out = bar_chart([Bar("a", None, "a: nothing")], "x", threshold=20, threshold_text="20x")
    assert ">NA<" in out and 'class="thresh"' in out


def _write_vcf(path, chroms=("chr1",), genomes=("a", "b", "c")):
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(genomes)
    records = [f"{chroms[0]}\t100\t.\tA\tG\t.\t.\t.\tGT\t1\t0\t0",
               f"{chroms[0]}\t2500\t.\tC\tT\t.\t.\t.\tGT\t1\t.\t1"]
    if len(chroms) > 1:
        records.append(f"{chroms[1]}\t50\t.\tG\tA,T\t.\t.\t.\tGT\t2\t1\t0")
    path.write_text("##fileformat=VCFv4.2\n" + header + "\n" + "\n".join(records) + "\n")


def test_read_vcf_counts_alt_and_missing_calls(tmp_path):
    _write_vcf(tmp_path / "s.vcf", ("c1", "c2"))
    snps, genomes = read_vcf(tmp_path / "s.vcf")
    assert genomes == 3 and [(s.chrom, s.pos, s.alt_count, s.missing) for s in snps] == \
        [("c1", 100, 1, 0), ("c1", 2500, 2, 1), ("c2", 50, 2, 0)]
    assert snps[2].ref == "G" and snps[2].alt == "A,T"


def test_genome_map_tracks_and_several_sequences(tmp_path):
    _write_vcf(tmp_path / "s.vcf", ("c1", "c<2>"))
    snps, genomes = read_vcf(tmp_path / "s.vcf")
    sequences = [("c1", 5000), ("c<2>", 1000)]
    out = genome_map(sequences, snps, genomes)
    assert out.count('class="snp-all"') == 3 and out.count('class="snp-miss"') == 2  # 2 + legend each
    assert "c1:100 A&gt;G: alternate allele in 1 of 3 genomes; called in every genome" in out
    assert "c1:2,500 C&gt;T: alternate allele in 2 of 3 genomes; 1 missing call" in out
    assert "c&lt;2&gt;" in out and "c<2>" not in out and "5,000 bp" in out and "1,000 bp" in out
    assert "nbar" not in out and "N bases" not in out
    tracks = {"c1": [0, 3, 0, 0, 0, 0], "c<2>": [12, 0]}
    out = genome_map(sequences, snps, genomes, tracks, 1000)
    assert out.count('class="nbar"') == 2 and "max 12 N per 1 kb" in out and "per 1 kb, all assemblies" in out
    assert "c1:1,001–2,000: 3 N bases" in out and "c&lt;2&gt;:1–1,000: 12 N bases" in out


def test_genome_map_bins_dense_snps():
    from bacon.report import Snp
    snps = [Snp("c", p, "A", "G", 1, 1 if p % 2 else 0) for p in range(1, 2001)]
    out = genome_map([("c", 2_000_000)], snps, 4)
    assert out.count("<line") < 100 and " SNPs at c:1–" in out and " with a missing call. 1 A&gt;G (1 alt, 1 missing)" in out
    assert out.count('class="snp-miss"') >= 2


def test_n_per_bin_matches_records_to_reference_sequences(tmp_path):
    (tmp_path / "s1.fasta").write_text(">s1_c1\nACGTNNNNAC\n>s1_c2 circular=true\nNNAC\n>s1_other\nNNNN\n")
    (tmp_path / "s2.fasta").write_text(">s2_c1\n" + "A" * 4 + "N" * 6 + "\n")
    counts = n_per_bin([tmp_path / "s1.fasta", tmp_path / "s2.fasta"], [("c1", 10), ("c2", 4)], 4)
    assert counts == {"c1": [0, 8, 2], "c2": [2, 0]}  # s1_other is not a reference sequence
    (tmp_path / "my_sample_2.fasta").write_text(">my_sample_2_c1\nNNNNACGTAC\n")  # "_" in the sample name
    counts = n_per_bin([tmp_path / "my_sample_2.fasta"], [("c1", 10), ("c2", 4)], 4)
    assert counts == {"c1": [4, 0, 0], "c2": [0, 0]}


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
    assert "in tree order" in hm["description"]
    assert distance_heatmap(dist)["xcats"] == ["Reference", "a", "b"]
    assert "the order of the distance table" in distance_heatmap(dist)["description"]
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


def test_report_escapes_the_reference_name_and_names_the_parsnp_tree(tmp_path):
    out = _folder(tmp_path, {"method": "parsnp", "core_snps": 3, "tree_method": "iqtree", "tree_svg": "t.svg",
                             "distances": "4_compared/ska/snp_distances.tsv"},
                  reference={"file": "/x/ref<1>&b.fa", "length": 100})
    (out / "t.svg").write_text("<svg></svg>")
    page = build_report(out)
    assert "ref&lt;1&gt;&amp;b.fa" in page and "ref<1>" not in page
    assert "Parsnp core-genome alignment; IQ-TREE, midpoint-rooted" in page
    assert '<meta name="color-scheme" content="light dark">' in page


def test_report_command_on_a_folder_that_is_not_bacon(tmp_path):
    import subprocess
    import sys
    done = subprocess.run([sys.executable, "-m", "bacon.report", str(tmp_path)], capture_output=True, text=True)
    assert done.returncode == 1 and "not a BACoN output folder" in done.stderr and "Traceback" not in done.stderr
    (tmp_path / "bacon.log").write_text("2026-10-07 [INFO] BACoN\n")  # A run that did not finish (no run_info)
    done = subprocess.run([sys.executable, "-m", "bacon.report", str(tmp_path)], capture_output=True, text=True)
    assert done.returncode == 1 and "the run did not finish: resume it to get a report" in done.stderr


TREE = "((a:0.001,Reference:0.0005)0.9:0.002,b:0.003);\n"


def _full_run(tmp_path, assembler="samtools", tree=TREE, vcf=True, reference=True, assemblies=True):
    comparison = {"method": "ska", "tree_method": "fasttree", "core_snps": 3, "tree": "4_compared/ska/tree.nwk",
                  "tree_svg": "4_compared/ska/tree.svg", "distances": "4_compared/ska/snp_distances.tsv",
                  "vcf": "4_compared/ska/snps.vcf"}
    summary = ("Sample\tStatus\tEst_depth\tN_bases\tAssembly_length\tNote\na\tok\t15.0\t2\t1000\tlow depth (15x)\n"
               "b\tok\t80.0\t0\t1000\t\n")
    matrix = "snp-dists\tReference\ta\tb\nReference\t0\t0\t3\na\t0\t0\t3\nb\t3\t3\t0\n"
    out = _folder(tmp_path, comparison, summary=summary, matrix=matrix,
                  reference={"file": "/x/ref.fa", "length": 1000, "sequences": 1})
    info = json.loads((out / "run_info.json").read_text())
    info["settings"] = {**SETTINGS, "assembler": assembler}
    (out / "run_info.json").write_text(json.dumps(info))
    compared = out / "4_compared" / "ska"
    if tree is not None:
        (compared / "tree.nwk").write_text(tree)
    (compared / "tree.svg").write_text("<svg><text>fallback tree</text></svg>")
    if vcf:
        _write_vcf(compared / "snps.vcf", ("c1",), ("a", "b"))
    if reference:
        (out / "reference.fasta").write_text(">c1 some description\n" + "ACGT" * 1000 + "\n")
    if assemblies:
        folder = out / "3_assembled" / "all_assemblies"
        folder.mkdir(parents=True)
        (folder / "a.fasta").write_text(">a_c1\n" + "N" * 10 + "ACGT" * 997 + "\n")
        (folder / "b.fasta").write_text(">b_c1\n" + "ACGT" * 1000 + "\n")
    return out


def test_report_draws_every_figure_for_a_templated_run(tmp_path):
    page = build_report(_full_run(tmp_path))
    for n in range(1, 6):
        assert f"<b>Figure {n}.</b>" in page
    assert "fallback tree" not in page and ">0.9</text>" in page and "ladderized" in page
    assert 'class="nbar"' in page and "templated assemblies per 1 kb" in page and "c1:1–1,000: 10 N bases" in page
    assert "c1:100 A&gt;G" in page and 'class="snp-miss"' in page
    assert ">15x</text>" in page and "20x flag" in page
    assert "group 1 (2)" in page and "a, Reference" in page  # a and Reference identical: the heatmap block
    assert '<meta name="color-scheme" content="light dark">' in page and "prefers-color-scheme:dark" in page
    assert "<script src" not in page and "http" not in page.split("</style>")[0]  # Self-contained


def test_report_omits_the_n_track_for_de_novo_assemblies(tmp_path):
    page = build_report(_full_run(tmp_path, assembler="flye"))
    assert 'class="nbar"' not in page and "No N track: flye assemblies are de novo" in page
    assert "<b>Figure 5.</b>" in page  # The SNP track is still drawn


def test_report_falls_back_to_the_drawn_tree_and_survives_missing_files(tmp_path):
    page = build_report(_full_run(tmp_path, tree="not a tree at all"))
    assert "fallback tree" in page and "midpoint-rooted; internal labels are supports" in page
    assert "in the order of the distance table" in page
    page = build_report(_full_run(tmp_path / "x", tree=None, vcf=False, reference=False, assemblies=False))
    assert "fallback tree" in page and "Genome map" not in page and "snps.vcf was not found" in page
    assert "<b>Figure 3.</b>" in page and "<b>Figure 5.</b>" not in page
    out = _full_run(tmp_path / "y", assemblies=False)
    page = build_report(out)
    assert "No N track: the assemblies were not found" in page and "<b>Figure 5.</b>" in page
    (out / "4_compared" / "ska" / "snps.vcf").write_text("garbage\n")
    assert "Genome map" in build_report(out)  # No records: an empty track, not an error
    (out / "reference.fasta").write_text("")
    assert "Genome map" not in build_report(out)


def test_report_escapes_names_in_every_figure(tmp_path):
    out = _full_run(tmp_path, tree="(('<x1>':0.001,Reference:0.0005)0.9:0.002,'<x2>':0.003);\n")
    (out / "summary.tsv").write_text("Sample\tStatus\tEst_depth\tN_bases\tAssembly_length\tNote\n<x1>\tok\t15.0\t2"
                                     "\t1000\tnote <n>\n<x2>\tok\t80.0\t0\t1000\t\n")
    matrix = "snp-dists\tReference\t<x1>\t<x2>\nReference\t0\t0\t3\n<x1>\t0\t0\t3\n<x2>\t3\t3\t0\n"
    (out / "4_compared" / "ska" / "snp_distances.tsv").write_text(matrix)
    _write_vcf(out / "4_compared" / "ska" / "snps.vcf", ("c<1>",), ("<x1>", "<x2>"))
    (out / "reference.fasta").write_text(">c<1>\n" + "ACGT" * 1000 + "\n")
    (out / "3_assembled" / "all_assemblies" / "<x1>.fasta").write_text("><x1>_c<1>\n" + "N" * 10 + "ACGT" * 997 + "\n")
    page = build_report(out)
    assert "<x1>" not in page and "<x2>" not in page and "<n>" not in page and "c<1>" not in page
    assert page.count("&lt;x1&gt;") > 8 and "note &lt;n&gt;" in page and "c&lt;1&gt;:1–1,000: 10 N bases" in page
    assert "<title>&lt;x1&gt; – &lt;x2&gt;: 3 SNPs</title>" in page and "group 1 (2)" in page


def test_genome_map_without_any_n_base(tmp_path):
    from bacon.report import _Figures, _genome_map_section
    out = tmp_path / "out"
    (out / "3_assembled" / "all_assemblies").mkdir(parents=True)
    (out / "3_assembled" / "all_assemblies" / "a.fasta").write_text(">a_chr\nACGTACGTAC\n")
    reference = out / "reference.fasta"
    reference.write_text(">chr\nACGTACGTAC\n")
    vcf = out / "snps.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\n"
                   "chr\t3\t.\tG\tT\t.\t.\t.\tGT\t1\n")
    html_text, _ = _genome_map_section(out, vcf, reference, {"assembler": "samtools"}, _Figures(), 1)
    assert "none of the 1 templated assemblies has an N base" in html_text and "N per" not in html_text


# ---------------------------------------------------------------------------------------------------------------
# Annotations: gene track, region band, SNP table
# ---------------------------------------------------------------------------------------------------------------

def _genbank(name="c1", length=4000, genes=True, repeats=True, gene_name="genA"):
    """A GenBank record of `length` bp of ACGT repeats: a + strand CDS at 1..300 (codon 2 is TAC: position 5 is
    its middle base), a - strand tRNA at 2000..2072, and two inverted repeats making LSC/IRb/SSC/IRa."""
    seq = ("ACGT" * (length // 4 + 1))[:length]
    features = ""
    if genes:
        features += (f"     gene            1..300\n                     /gene=\"{gene_name}\"\n"
                     f"     CDS             1..300\n                     /gene=\"{gene_name}\"\n"
                     "                     /codon_start=1\n                     /transl_table=11\n"
                     "                     /product=\"protein A\"\n"
                     "     gene            complement(2000..2072)\n                     /gene=\"trnQ\"\n"
                     "     tRNA            complement(2000..2072)\n                     /gene=\"trnQ\"\n")
    if repeats:
        features += ("     repeat_region   1001..1600\n                     /rpt_type=inverted\n"
                     "     repeat_region   2401..3000\n                     /rpt_type=inverted\n")
    origin = "\n".join(f"{i + 1:>9} " + " ".join(seq[j:j + 10].lower() for j in range(i, min(i + 60, length), 10))
                       for i in range(0, length, 60))
    return (f"LOCUS       {name}  {length} bp    DNA     circular PLN 01-JAN-2026\nDEFINITION  A test.\n"
            f"VERSION     {name}\nFEATURES             Location/Qualifiers\n{features}ORIGIN      \n{origin}\n//\n")


def _annotated_run(tmp_path, **kw):
    out = _full_run(tmp_path)
    (out / "4_compared" / "ska" / "snps.vcf").write_text(  # 101 is the middle base of codon 34 (TAC) of the CDS
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n"
        "c1\t101\t.\tA\tG\t.\t.\t.\tGT\t1\t0\nc1\t2500\t.\tT\tC\t.\t.\t.\tGT\t1\t.\n")
    (out / "annotation.gb").write_text(_genbank(**kw))
    info = json.loads((out / "run_info.json").read_text())
    info["annotation"] = {"file": "/x/NC_1.gb", "format": "genbank", "copy": "annotation.gb", "transl_tables": [11]}
    (out / "run_info.json").write_text(json.dumps(info))
    return out


def test_report_with_annotation_draws_genes_regions_and_the_snp_table(tmp_path):
    page = build_report(_annotated_run(tmp_path))
    assert 'class="gene-cds"' in page and 'class="gene-rna"' in page and 'class="band-ir"' in page
    assert "LSC 2,000 bp" in page and "IRb 600 bp" in page and "SSC 800 bp" in page and "IRa 600 bp" in page
    assert "<b>Figure 5.</b>" in page and "<b>Figure 6.</b>" not in page  # The map is still one figure
    assert "Genes from <b>NC_1.gb</b> (2)" in page  # The file as given, not BACoN's copy
    assert "The band above the genes shows the LSC/IRb/SSC/IRa regions from the annotated inverted repeats." in page
    assert "The LSC/IRb/SSC/IRa regions were derived from the annotated inverted repeats." in page  # Methods
    # The SNP at 101 (A>G) is in the CDS: codon 34 is TAC (positions 100..102), so A>G gives TGC: Y34C
    assert "<h3>SNPs</h3>" in page and 'class="snps sortable"' in page
    assert '<td class="gene" data-v="genA">genA</td>' in page and ">TAC&gt;TGC<" in page and ">Y34C<" in page
    assert ">missense<" in page and "2 SNPs on the annotated sequences: 1 in IRa, 1 in the LSC" in page
    assert "1 in coding sequences (1 missense); 1 intergenic" in page
    assert "intergenic between trnQ and genA" in page  # 2,500 is in IRa, after the tRNA, wrapping to genA
    assert "genA (CDS) · missense Y34C (TAC&gt;TGC)" in page  # The SNP tick's hover
    assert "read from the annotation NC_1.gb" in page and "translation table 11" in page  # Methods
    assert 'class="t-tiny t-ink2 t-gene"' not in page  # genA has one SNP: no label


def test_report_without_annotation_has_no_gene_track(tmp_path):
    page = build_report(_full_run(tmp_path))
    assert 'class="gene-' not in page and "<h3>SNPs</h3>" not in page and "Genes from" not in page


def test_report_escapes_gene_names_everywhere(tmp_path):
    out = _annotated_run(tmp_path, gene_name="g<1>&")
    (out / "4_compared" / "ska" / "snps.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n"
        "c1\t5\t.\tA\tC\t.\t.\t.\tGT\t1\t0\nc1\t9\t.\tA\tC\t.\t.\t.\tGT\t1\t0\n")
    page = build_report(out)
    assert "g<1>&" not in page and page.count("g&lt;1&gt;&amp;") >= 5  # Gene rect, label, table, hover, summary
    assert 'class="t-tiny t-ink2 t-gene">g&lt;1&gt;&amp;' in page  # Two SNPs: labelled


def test_report_with_a_broken_annotation_keeps_the_map(tmp_path):
    out = _annotated_run(tmp_path)
    (out / "annotation.gb").write_text("garbage\n")
    page = build_report(out)
    assert "Genome map" in page and "could not be used" in page and 'class="gene-' not in page
    (out / "annotation.gb").unlink()  # The copy is gone (an old folder): no annotation, no error
    page = build_report(out)
    assert "Genome map" in page and "could not be used" not in page


def test_report_finds_the_annotation_copy_without_run_info(tmp_path):
    from bacon.report import annotation_file
    out = _annotated_run(tmp_path)
    assert annotation_file(out, {}) == out / "annotation.gb"
    assert annotation_file(out, {"annotation": {"copy": "annotation.gff3"}}) == out / "annotation.gb"
    assert annotation_file(tmp_path, {}) is None


def test_annotation_warnings_are_shown(tmp_path):
    out = _annotated_run(tmp_path, name="other", length=4001)
    page = build_report(out)
    assert "match no reference sequence by name" in page and 'class="gene-' not in page


def test_genome_map_merges_genes_when_there_are_too_many(tmp_path, monkeypatch):
    import bacon.report as report
    from bacon.annotation import load_annotation
    gb = tmp_path / "a.gb"
    gb.write_text(_genbank())
    annotation = load_annotation(gb, [("c1", 4000)])
    monkeypatch.setattr(report, "MAX_GENE_RECTS", 1)
    out = genome_map([("c1", 4000)], [], 2, annotation=annotation)
    assert 'class="gene-dense"' in out and "gene-cds" not in out and "t-gene" not in out
    assert "1 gene: genA" in out and "1 gene: trnQ" in out
    page = build_report(_annotated_run(tmp_path))
    assert "merged per pixel" in page


def test_place_labels_alternates_rows_and_caps():
    from bacon.annotation import Gene
    from bacon.report import MAX_GENE_LABELS, place_labels
    genes = [Gene("c", f"gene{i:03d}", "CDS", 1, [(i * 10 + 1, i * 10 + 10)], snps=2) for i in range(100)]
    genes[50].snps = 1  # Not labelled
    placed, rows = place_labels(genes, 100, 1100, 1.0)  # 10 px per gene: labels overlap on one row
    assert rows == 2 and len(placed) <= MAX_GENE_LABELS and all(g.snps >= 2 for _, _, g in placed)
    assert [row for _, row, _ in placed][:2] == [0, 1]
    assert place_labels(genes, 100, 1100, 0.01)[1] <= 2
    assert place_labels([], 0, 100, 1.0) == ([], 0)


def test_snp_table_caps_rows(tmp_path, monkeypatch):
    import bacon.report as report
    from bacon.report import Snp, snp_table
    monkeypatch.setattr(report, "MAX_SNP_ROWS", 2)
    snps = [Snp("c1", p, "A", "G", 1, 0) for p in (1, 2, 3)]
    out = snp_table(snps, {}, 2)
    assert out.count("<tr>") == 3 + 1 - 1 and "The first 2 of the 3 SNPs" in out  # Header row + 2 rows
    assert out.count("<th ") == 10 and snp_table(snps, {}, 2, regions=False).count("<th ") == 9  # No Region column


def test_snp_summary_says_when_overlapping_genes_count_a_snp_twice():
    from types import SimpleNamespace

    from bacon.report import Snp, snp_summary
    snps = [Snp("c", 5, "A", "G", 1, 0), Snp("c", 9, "A", "G", 1, 0)]
    info = {("c", 5): SimpleNamespace(region="", context="pseudogene / CDS", effects=[], genes=[]),
            ("c", 9): SimpleNamespace(region="", context="CDS", effects=[], genes=[])}
    text = snp_summary(snps, info, SimpleNamespace(sequences={}))
    assert "2 in coding sequences" in text and "counted for each" in text
    info[("c", 5)].context = "CDS"
    assert "counted for each" not in snp_summary(snps, info, SimpleNamespace(sequences={}))


# ---------------------------------------------------------------------------------------------------------------
# Sample metadata in the report
# ---------------------------------------------------------------------------------------------------------------

def _colouring(values: dict[str, str], column="group"):
    from bacon.metadata import Metadata
    from bacon.report import colouring
    return colouring(Metadata([column], {k: {column: v} for k, v in values.items()}), column)


def test_colouring_orders_values_and_colours_them():
    c = _colouring({"a": "10", "b": "9", "c": "", "d": "9"})
    assert c.values == ["9", "10"] and c.slot("b") == 0 and c.slot("a") == 1 and c.slot("c") is None
    assert c.colour("b") == "var(--s1)" and c.colour("c") == "var(--muted)" and c.colour("zz") == "var(--muted)"
    assert c.title("a") == "a: group 10" and c.title("c") == "c: no group"


def test_heatmap_second_band_and_legend():
    names = ["a", "b", "c", "<d>"]
    m = {x: {y: 0 if x == y else 12 for y in names} for x in names}
    m["a"]["b"] = m["b"]["a"] = 0
    groups = identical_groups(names, m)
    colours = _colouring({"a": "x<1", "b": "y", "c": "x<1", "<d>": ""}, column="Site<s>")
    out = heatmap(names, m, groups, colours)
    assert out.count("<title>a: Site&lt;s&gt; x&lt;1</title>") == 2 and out.count("<title>a: group 1</title>") == 2
    assert out.count("<title>&lt;d&gt;: no Site&lt;s&gt;</title>") == 2 and 'fill="var(--muted)"' in out
    assert "Site&lt;s&gt;:</text>" in out and "x&lt;1 (2)</text>" in out and "y (1)</text>" in out
    assert "no value (1)</text>" in out and "<d>" not in out and "x<1" not in out
    plain = heatmap(names, m, groups)
    assert "Site" not in plain and "no value" not in plain and len(plain) < len(out)
    # The value band is outside the group band: at the far left, the group band 14 px from the cells
    left_band = min(int(x) for x in __import__("re").findall(r'<rect x="(\d+)" y="\d+" width="10"', out))
    assert f'<rect x="{left_band}"' in out and f'<rect x="{left_band + 14}"' in out


def test_tree_svg_with_values():
    root = parse("((a:0.01,b:0.01)0.95:0.02,(Reference:0.005,'<c>':0.005)0.80:0.03);")
    colours = _colouring({"a": "x", "b": "y", "<c>": "<v>"})
    out = tree_svg(root, {"a": 0, "b": 0}, "ref.fa", 100, colours)
    assert out.count("<circle") == 4 and 'r="4" fill="var(--s2)"><title>a: group x</title>' in out  # <v>, x, y
    assert 'fill="none" stroke="var(--muted)"><title>Reference: no group</title>' in out  # Hollow: no value
    assert '>a <tspan class="t-muted" font-weight="400">x</tspan></text>' in out
    assert '>&lt;c&gt; <tspan class="t-muted" font-weight="400">&lt;v&gt;</tspan>' in out and "<v>" not in out
    assert "Reference <tspan" in out and "ref.fa" in out
    plain = tree_svg(root, {"a": 0, "b": 0}, "ref.fa", 100)
    assert "<circle" not in plain and "t-muted\" font-weight=\"400\">x<" not in plain


def test_cross_table_counts():
    from bacon.report import cross_table
    names = ["Reference", "a", "b", "c", "d", "<e>"]
    groups = [["Reference", "a", "b"], ["c", "d"]]
    colours = _colouring({"a": "x", "b": "y", "c": "x", "d": "x", "<e>": "<y>"}, column="g")
    out = cross_table(groups, names, colours)
    rows = out.split("<tr>")[1:]
    assert "&lt;y&gt;</th>" in rows[0] and "no value</th>" in rows[0] and "Total</th>" in rows[0]
    assert rows[1].count("<td") == 6 and ">1</td>" in rows[1] and ">3</td>" in rows[1]  # group 1: x 1, y 1, none 1
    assert ">2</td>" in rows[2] and rows[2].count('class="num zero">–</td>') == 3  # group 2: x 2
    assert "not in a group" in rows[3] and ">1</td>" in rows[3] and "<e>" not in out
    without_missing = cross_table([["a", "b"]], ["a", "b"], colours)
    assert "no value" not in without_missing and "not in a group" not in without_missing


def _metadata_run(tmp_path, copy, color_by="group", record=True):
    out = _full_run(tmp_path)
    (out / "metadata.tsv").write_text(copy)
    if record:
        info = json.loads((out / "run_info.json").read_text())
        info["metadata"] = {"file": "/x/meta <m>.tsv", "copy": "metadata.tsv", "columns": copy.split("\n")[0].split("\t")[1:],
                            "sample_sheet_columns": [], "color_by": color_by}
        (out / "run_info.json").write_text(json.dumps(info))
    return out


def test_report_with_metadata_table_figures_and_escaping(tmp_path):
    out = _metadata_run(tmp_path, "sample\tgroup\tn <x>\tnote\na\tG<1>\t10\tsome text\nb\t\t9.5\t\n")
    page = build_report(out)
    head = page.split("<table class=\"samples sortable\">")[1].split("</thead>")[0]
    assert head.index(">Sample<") < head.index('<th class="md">group</th>') < head.index('class="num md">n &lt;x&gt;<') \
        < head.index(">Status<")
    assert '<td class="num md" data-v="10">10</td>' in page and '<td class="num md" data-v="9.5">9.5</td>' in page
    assert '<td class="md" data-v="G&lt;1&gt;"><span class="swatch" style="background:var(--s1)"></span>G&lt;1&gt;' \
        in page
    assert '<td class="md" data-v=""></td>' in page and "G<1>" not in page and "n <x>" not in page
    assert "Metadata from <b>meta &lt;m&gt;.tsv</b>: 3 columns (group, n &lt;x&gt;, note); 2 of 2 samples" in page
    assert "coloured by <b>group</b>" in page and "give each genome's <b>group</b>" in page
    assert "<title>a: group G&lt;1&gt;</title>" in page and "G&lt;1&gt; (1)</text>" in page  # Band and legend
    assert "no value (2)</text>" in page  # b and the Reference
    assert "Genomes of each group by <b>group</b>" in page and '<table class="cross">' in page
    assert ">Metadata</dt><dd><code>/x/meta &lt;m&gt;.tsv (copy metadata.tsv; colours by group)</code>" in page


def test_report_metadata_without_colours_and_fallbacks(tmp_path):
    out = _metadata_run(tmp_path, "sample\tgroup\na\tx\nb\ty\n", color_by=None)
    page = build_report(out)
    assert "No column colours the figures" in page and "<circle" not in page and "group:</text>" not in page
    assert '<table class="cross">' not in page and '<th class="md">group</th>' in page
    # No record in run_info (a hand-made copy): the first usable column
    out = _metadata_run(tmp_path / "x", "sample\tgroup\na\tx\nb\ty\n", record=False)
    page = build_report(out)
    assert "coloured by <b>group</b>" in page and "Metadata from <b>metadata.tsv</b> (found in the output folder; " \
        "not given to this run, so run_info.json and the MultiQC table do not have it): 1 column" in page
    # A recorded column missing from the copy: no colours; an unreadable copy: a note
    out = _metadata_run(tmp_path / "y", "sample\tother\na\tx\nb\ty\n", color_by="group")
    assert "No column colours the figures" in build_report(out)
    (out / "metadata.tsv").write_text("name\tother\na\tx\n")
    page = build_report(out)
    assert "No metadata: Metadata file" in page and "needs a &#x27;sample&#x27; column" in page
    assert '<th class="md">' not in page


def test_report_without_metadata_is_unchanged(tmp_path):
    out = _full_run(tmp_path)
    page = build_report(out)
    assert "Metadata" not in page and "<circle" not in page and 'class="md' not in page
    assert 'class="cross"' not in page and "no value" not in page
    info = json.loads((out / "run_info.json").read_text())
    info["metadata"] = None  # As a run without metadata records it
    (out / "run_info.json").write_text(json.dumps(info))
    assert build_report(out) == page


def test_multiqc_table_with_metadata():
    from bacon.metadata import Metadata
    metadata = Metadata(["Group", "Site name", "Status"], {"s1": {"Group": "A", "Site name": "n", "Status": "x"}})
    table = sample_table(ROWS, "run", metadata)
    keys = list(table["headers"])
    assert keys[:3] == ["meta_group", "meta_site_name", "meta_status"] and keys[3] == "Status"
    assert list(sample_table(ROWS, "run", Metadata(["a b", "a_b"], {}))["headers"])[:2] == ["meta_a_b", "meta_a_b_2"]
    assert table["headers"]["meta_site_name"] == {"title": "Site name", "description": "Metadata: Site name"}
    assert table["data"]["s1"]["meta_group"] == "A" and "meta_group" not in table["data"]["s2"]
    assert table["description"].endswith("with the sample metadata.")


def test_with_a_colour_column_colour_means_the_metadata_only():
    from bacon.report import Colouring, heatmap, tree_svg
    colours = Colouring("group", ["A", "B"], {"a": "A", "b": "B", "c": "A"})
    root = parse("((a:0.01,b:0.01)0.95:0.02,(Reference:0.005,c:0.005)0.80:0.03);")
    tree = tree_svg(root, {"a": 0, "b": 0}, "ref.fa", 100, colours)
    assert "<title>group 1</title>" not in tree  # No group square: the circles carry the colour
    names = ["Reference", "a", "b", "c"]
    matrix = {x: {y: (0 if {x, y} <= {"a", "b"} or x == y else 3) for y in names} for x in names}
    page = heatmap(names, matrix, [["a", "b"]], colours)
    assert 'fill="var(--ink2)"><title>a: group 1</title>' in page  # Grey group band
    assert 'fill="var(--s1)"><title>a: group A</title>' in page  # Coloured metadata band
    assert 'fill="var(--s1)"><title>a: group 1</title>' in heatmap(names, matrix, [["a", "b"]])  # Without metadata


def test_grey_group_bands_alternate_along_the_axes():
    from bacon.report import Colouring, heatmap
    names = ["a", "b", "c", "d"]
    matrix = {x: {y: (0 if (x in "ab") == (y in "ab") else 4) for y in names} for x in names}
    colours = Colouring("group", ["A"], {n: "A" for n in names})
    page = heatmap(names, matrix, [["c", "d"], ["a", "b"]], colours)  # Group 2 comes first on the axes
    assert 'fill="var(--ink2)"><title>a: group 2</title>' in page
    assert 'fill="var(--muted)"><title>c: group 1</title>' in page  # Adjacent groups differ


def test_grey_group_bands_alternate_per_run_along_the_axes():
    from bacon.report import Colouring, heatmap
    names = ["a", "b", "c", "d", "e", "f"]
    groups = [["a", "c"], ["b", "e"], ["d", "f"]]  # Every group is split by the order of the axes
    matrix = {x: {y: (0 if any({x, y} <= set(g) for g in groups) or x == y else 4) for y in names} for x in names}
    colours = Colouring("group", ["A"], {n: "A" for n in names})
    page = heatmap(names, matrix, groups, colours)
    fills = [re.search(rf'fill="(var\(--\w+\))"><title>{n}: group \d</title>', page).group(1) for n in names]
    assert fills == ["var(--ink2)", "var(--muted)"] * 3  # Adjacent runs never share a grey (d after c)
    blocks = re.findall(r'width="4" height="\d+" fill="(var\(--\w+\))"><title>', page)  # The blocks on the right
    assert blocks == fills


def test_cross_table_columns_only_for_values_of_the_genomes_shown():
    from bacon.report import cross_table
    colours = _colouring({"a": "A", "b": "B", "z": "Z"})  # z failed: it is in no figure
    out = cross_table([["a", "b"]], ["Reference", "a", "b"], colours)
    head = out.split("</thead>")[0]
    assert "Z</th>" not in head and "A</th>" in head and "B</th>" in head and "no value</th>" in head
    assert 'style="background:var(--s2)"></span>B</th>' in head  # B keeps its slot colour of the legend
    assert out.split("<tr>")[2].count("<td") == 5  # group 1: A, B, no value, Total


def test_metadata_columns_named_like_the_tables_own_are_marked(tmp_path):
    out = _metadata_run(tmp_path, "sample\tStatus\tnote\tsite\ta\tfine\thello\tnorth\nb\tbad\t\tsouth\n".replace(
        "site\ta", "site\na"), color_by="site")
    page = build_report(out)
    head = page.split("<table class=\"samples sortable\">")[1].split("</thead>")[0]
    assert re.findall(r"<th[^>]*>([^<]*)</th>", head)[:5] == ["Sample", "Status (metadata)", "note (metadata)",
                                                              "site", "Status"]
    assert "Metadata from <b>meta &lt;m&gt;.tsv</b>: 3 columns (Status, note, site)" in page  # The file's names
    from bacon.metadata import Metadata
    table = sample_table(ROWS, "run", Metadata(["Status", "Depth", "sample", "site"], {}))
    assert [h["title"] for h in list(table["headers"].values())[:4]] == \
        ["Status (metadata)", "Depth (metadata)", "sample (metadata)", "site"]


def test_print_stylesheet_unrolls_tall_tables_and_wraps_cells(tmp_path):
    assert '<div class="tablewrap tall"><table class="snps' in build_report(_annotated_run(tmp_path))
    from bacon.report import CSS
    print_css = CSS.split("@media print{", 1)[1]
    assert ".tablewrap.tall{max-height:none;overflow:visible}" in print_css  # Every SNP row is printed
    assert "th,td{white-space:normal" in print_css and "th{position:static}" in print_css  # Wide tables wrap


def test_tree_svg_keeps_the_support_of_a_root_child_inside():
    from bacon.newick import midpoint_root
    from bacon.report import _text_px
    out = tree_svg(midpoint_root(parse("(((A:1,B:1)0.9:1,C:1)0.8:1,D:1,E:1);")), {}, "", None)
    supports = re.findall(r'<text x="([-\d.]+)" y="[-\d.]+" text-anchor="end" class="t-tiny t-muted">([^<]+)<', out)
    assert {s for _, s in supports} == {"0.9", "0.8"}
    assert all(float(x) - _text_px(s, 9.5) >= 0 for x, s in supports)  # Not off the left edge
    assert float(re.search(r'<path d="M([\d.]+),', out).group(1)) >= 16


def _blocky_matrix():
    """150 genomes at cell size 10: ten groups of 12 whose members alternate (blocks of one row), an eleventh with
    eleven members in a row and one at the end, and 18 singletons."""
    names = [f"n{i:03d}" for i in range(150)]
    group = {names[p]: p % 10 for p in range(120)}
    group.update({names[p]: 10 for p in range(120, 131)})
    group[names[149]] = 10
    matrix = {a: {b: 0 if a == b or group.get(a, -1) == group.get(b, -2) else 5 for b in names} for a in names}
    return names, matrix


def test_heatmap_labels_only_the_tall_blocks_and_fits_the_longest_label():
    from bacon.report import _text_px
    names, matrix = _blocky_matrix()
    groups = identical_groups(names, matrix)
    assert len(groups) == 11 and groups[10][0] == "n120"
    out = heatmap(names, matrix, groups)
    labels = re.findall(r'<text x="(\d+)" y="\d+" class="t-small t-ink">([^<]+)</text>', out)
    assert [label for _, label in labels] == ["group 11 (11 of 12)"]  # The one block of 12 px or more
    assert out.count("<title>group ") == 122  # Every block has its bar, with the label on hover
    width = float(re.search(r'viewBox="0 0 (\d+)', out).group(1))
    assert width - float(labels[0][0]) >= _text_px(labels[0][1], 11)


def test_bar_chart_cuts_long_names_and_keeps_them_on_hover():
    name = "sample_" + "x" * 53  # 60 characters
    out = bar_chart([Bar(name, 5, "h"), Bar("short", 3, "h")], "t")
    assert f">{name}</text>" not in out and f"…<title>{name}</title></text>" in out and ">short</text>" in out
    label = re.search(r'text-anchor="end" class="t-small t-ink2">([^<]*)…', out).group(1)
    assert name.startswith(label) and 0.58 * 11 * (len(label) + 1) <= 220 - 16


def test_bar_chart_of_zeros_has_no_made_up_scale():
    rows = [{"Sample": "a", "Status": "ok", "N_bases": "0"}, {"Sample": "b", "Status": "ok", "N_bases": "0"}]
    out = bar_chart(n_bars(rows), "N bases", integers=True)
    ticks = re.findall(r'text-anchor="middle" class="t-tiny t-muted">([^<]*)</text>', out)
    assert ticks == ["0"] and out.count('class="bar0"') == 2


def test_n_per_bin_rescales_an_assembly_of_another_length(tmp_path):
    # A 100-bp insertion before a run of N shifts it in the consensus; a 100-bp deletion pulls it forward
    (tmp_path / "ins.fasta").write_text(">ins_chr1\n" + "A" * 2100 + "N" * 100 + "A" * 400 + "\n")  # 2,600 bp
    (tmp_path / "del.fasta").write_text(">del_chr1\n" + "A" * 1000 + "N" * 100 + "A" * 1300 + "\n")  # 2,400 bp
    (tmp_path / "edge.fasta").write_text(">edge_chr1\n" + "A" * 450 + "N" * 100 + "A" * 450 + "\n")  # 1,000 bp
    assert n_per_bin([tmp_path / "ins.fasta"], [("chr1", 2500)], 500) == {"chr1": [0, 0, 0, 0, 100, 0]}
    assert n_per_bin([tmp_path / "del.fasta"], [("chr1", 2500)], 500) == {"chr1": [0, 0, 100, 0, 0, 0]}
    # 2,101–2,200 x 2500/2600 = 2,020–2,115 (the fifth bin); 1,001–1,100 x 2500/2400 = 1,042–1,146 (the third)
    assert n_per_bin([tmp_path / "edge.fasta"], [("chr1", 1040)], 520) == {"chr1": [50, 50, 0]}  # Split, same sum
    both = n_per_bin([tmp_path / "ins.fasta", tmp_path / "del.fasta"], [("chr1", 2500)], 500)
    assert sum(both["chr1"]) == 200
    # Rescaling is for indels: up to 5% of the length. A record of another length is a partial consensus, counted
    # at its own positions (rescaled, the half-length record's N run would land in the last bin)
    (tmp_path / "near.fasta").write_text(">near_chr1\n" + "A" * 890 + "N" * 10 + "A" * 60 + "\n")  # 960 bp: 4%
    assert n_per_bin([tmp_path / "near.fasta"], [("chr1", 1000)], 100)["chr1"][8:10] == [0, 10]  # 891–900 x 1.04
    (tmp_path / "half.fasta").write_text(">half_chr1\n" + "A" * 490 + "N" * 10 + "\n")  # 500 bp of 1,000
    assert n_per_bin([tmp_path / "half.fasta"], [("chr1", 1000)], 100) == {"chr1": [0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0]}


def test_read_vcf_classes_every_genotype_form(tmp_path):
    forms = ["0/1", "1|0", "./.", "0:35", ".", "2", "0", "1", "1:12:3", ".:0", "0/0", "", "10", "0|.", "./1", "1/."]
    genomes = [f"g{i}" for i in range(len(forms))]
    (tmp_path / "m.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(genomes)
        + "\nc\t5\t.\tA\tG,T\t.\t.\t.\tGT:DP\t" + "\t".join(forms) + "\nc\t9\t.\tA\tG\t.\t.\t.\tGT\t"
        + "\t".join("0" for _ in forms) + "\n")
    snps, count = read_vcf(tmp_path / "m.vcf")
    assert count == 16 and [(s.alt_count, s.missing) for s in snps] == [(8, 4), (0, 0)]

    def slow(gt):  # The one-regex-per-call classification
        alleles = re.split(r"[/|]", gt.split(":")[0])
        return ("missing" if all(x in (".", "") for x in alleles)
                else "alt" if any(x.isdigit() and int(x) > 0 for x in alleles) else "ref")

    from bacon.report import _call_kind
    assert [_call_kind(gt) for gt in forms] == [slow(gt) for gt in forms]


def test_report_tile_has_a_thousands_separator(tmp_path):
    out = _folder(tmp_path, {"method": "ska", "core_snps": 12345, "distances": "4_compared/ska/snp_distances.tsv"})
    assert '<div class="tile-v">12,345</div>' in build_report(out)


def test_genome_map_axis_in_mb_for_long_sequences():
    out = genome_map([("c", 5_000_000)], [], 2)
    ticks = re.findall(r'text-anchor="middle" class="t-tiny t-muted">([^<]*)</text>', out)
    assert ticks[:2] == ["0 Mb", "0.5 Mb"] and ticks[-1] == "5 Mb"
    assert "20 kb" in genome_map([("c", 200_000)], [], 2)


def _wcag(a, b):
    def lum(h):
        h = h if len(h) == 7 else "#" + "".join(c * 2 for c in h[1:])
        r, g, b = (int(h[i:i + 2], 16) / 255 for i in (1, 3, 5))
        f = (lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4)
        return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)
    la, lb = lum(a), lum(b)
    return (max(la, lb) + 0.05) / (min(la, lb) + 0.05)


def test_small_text_and_in_cell_text_reach_wcag_contrast():
    from bacon.report import MISSING_COLOUR, TOKENS_DARK, TOKENS_LIGHT
    missing = re.fullmatch(r"var\(--(\w+)\)", MISSING_COLOUR).group(1)
    for tokens in (TOKENS_LIGHT, TOKENS_DARK):
        t = dict(re.findall(r"--([\w-]+):(#[0-9a-fA-F]{3,6})", tokens))
        assert _wcag(t["muted"], t["surface"]) >= 4.5 and _wcag(t["muted"], t["bg"]) >= 4.5  # Small muted text
        assert all(_wcag(t[f"qt{q}"], t[f"q{q}"]) >= 4.5 for q in range(5))  # Distances written in the cells
        assert _wcag(t[missing], t["surface"]) >= 3  # The "no value" band, hollow circle and swatch


def test_gene_title_map_and_labels_use_the_ranges_of_a_gene_in_pieces():
    from bacon.annotation import Gene, SequenceAnnotation
    from bacon.report import _gene_rows, _gene_title, place_labels
    wrapped = Gene("c", "psbA", "CDS", 1, [(1, 33), (1990, 2000)], snps=2)  # join(1990..2000,1..33), 2,000 bp
    spliced = Gene("c", "rps12", "CDS", -1, [(100, 130), (1500, 1600)], snps=2)  # Trans-spliced
    assert "1–33 + 1,990–2,000 (44 bp)" in _gene_title(wrapped) and "1–2,000" not in _gene_title(wrapped)
    assert "100–130 + 1,500–1,600 (132 bp)" in _gene_title(spliced)
    seq_ann = SequenceAnnotation("c", 2000, [wrapped, spliced], [], True)
    dense = "".join(_gene_rows(seq_ann, 0, 2000, 1.0, 0, True))
    boxes = [(int(x), int(w)) for x, w in re.findall(r'<rect x="(\d+)\.0" y="\d+" width="(\d+)"', dense)]
    assert all(w <= 110 for _, w in boxes) and len(boxes) == 4  # The pieces, not the hull of each gene
    placed, _ = place_labels([wrapped, spliced], 0, 2000, 1.0)
    centres = {g.name: x for x, _, g in placed}
    assert centres["psbA"] <= 40 and 1540 <= centres["rps12"] <= 1560  # On the largest piece, not the hull's middle


def test_heatmap_classes_at_the_upper_edge_of_each_bin():
    import re

    from bacon.report import distance_bins, heatmap
    bins = distance_bins(118)
    assert bins == [(0, 0), (1, 3), (4, 10), (11, 30), (31, 118)]
    names = ["Reference", "a", "b", "c", "d"]
    pairs = {("a", "b"): 3, ("a", "c"): 10, ("a", "d"): 30, ("b", "c"): 118}
    matrix = {x: {y: 0 for y in names} for x in names}
    for (x, y), d in pairs.items():
        matrix[x][y] = matrix[y][x] = d
    page = heatmap(names, matrix)
    groups: dict[str, str] = {}
    for cls, content in re.findall(r'<g class="(q\d)"[^>]*>(.*?)</g>', page, re.S):  # Cells (and legend swatches)
        groups[cls] = groups.get(cls, "") + content
    for (x, y), d, cls in [(("a", "b"), 3, "q1"), (("a", "c"), 10, "q2"), (("a", "d"), 30, "q3")]:
        assert f"{x} – {y}: {d} SNPs" in groups[cls], (d, cls)


def test_snp_table_names_the_sequence_when_the_reference_has_several(tmp_path):
    out = _annotated_run(tmp_path)
    page = build_report(out)
    assert '<th class="">Sequence</th>' not in page and "on sequences without annotation" not in page
    (out / "reference.fasta").write_text(">c1 some description\n" + "ACGT" * 1000 + "\n>c2\n" + "ACGT" * 100 + "\n")
    (out / "4_compared" / "ska" / "snps.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n"
        "c1\t101\t.\tA\tG\t.\t.\t.\tGT\t1\t0\nc1\t2500\t.\tT\tC\t.\t.\t.\tGT\t1\t.\n"
        "c2\t50\t.\tG\tA\t.\t.\t.\tGT\t1\t0\n")
    page = build_report(out)
    assert '<th class="">Sequence</th><th class="num">Position</th>' in page
    assert '<tr><td class="" data-v="c2">c2</td><td class="num" data-v="50">50</td>' in page
    assert "2 SNPs on the annotated sequences: 1 in IRa, 1 in the LSC. 1 SNP is on sequences without annotation " \
        "(no gene, context or effect). 1 in coding sequences" in page


def test_report_of_a_copied_folder_reads_its_own_files(tmp_path):
    import shutil
    out = _full_run(tmp_path)
    info = json.loads((out / "run_info.json").read_text())
    info["settings"]["output"] = str(out)  # As a run records them: absolute paths in the run's folder
    info["comparison"] = {k: str(out / v) if isinstance(v, str) and v.startswith("4_compared") else v
                          for k, v in info["comparison"].items()}
    (out / "run_info.json").write_text(json.dumps(info))
    copy = shutil.copytree(out, tmp_path / "copy")
    (copy / "4_compared" / "ska" / "snps.vcf").unlink()
    page = build_report(copy)
    assert "snps.vcf was not found" in page and "Genome map" not in page  # Not the original folder's file
    assert "SNP distances" in page and "not found" not in page.split("snps.vcf was not found")[0]
    assert "Genome map" in build_report(out)


def test_read_vcf_classes_each_genotype_once_whatever_the_other_fields(tmp_path, monkeypatch):
    import bacon.report as report_module
    classed = []
    real = report_module._call_kind
    monkeypatch.setattr(report_module, "_call_kind", lambda gt: classed.append(gt) or real(gt))
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n"
    rows = "".join(f"c1\t{pos}\t.\tA\tG\t.\t.\t.\tGT:DP\t1:{pos}\t0:{pos + 1}\t.:{pos + 2}\n" for pos in (10, 20, 30))
    (tmp_path / "s.vcf").write_text("##fileformat=VCFv4.2\n" + header + rows)
    snps, genomes = read_vcf(tmp_path / "s.vcf")
    assert genomes == 3 and [(s.alt_count, s.missing) for s in snps] == [(1, 1)] * 3
    assert sorted(classed) == [".", "0", "1"]  # Not once per distinct sample column (nine here)


# ---------------------------------------------------------------------------------------------------------------
# Regions detected in the sequence, and gene labels
# ---------------------------------------------------------------------------------------------------------------

def _plastid_genbank(seq: str, features: str, name: str = "c1") -> str:
    origin = "\n".join(f"{i + 1:>9} " + " ".join(seq[j:j + 10].lower() for j in range(i, min(i + 60, len(seq)), 10))
                       for i in range(0, len(seq), 60))
    return (f"LOCUS       {name}  {len(seq)} bp    DNA     circular PLN 01-JAN-2026\nDEFINITION  A plastid.\n"
            f"VERSION     {name}\nFEATURES             Location/Qualifiers\n{features}ORIGIN      \n{origin}\n//\n")


def _plastid_run(tmp_path, features: str | None = None, positions=(100, 12000), **kw):
    """A run on a 25 kb plastome-like reference (LSC 1-9,000, IRb 9,001-15,000, SSC, IRa 19,001-25,000) with SNPs
    at `positions`, and a GenBank annotation with `features` when given (None: no annotation)."""
    from tests.test_annotation import make_plastome
    out = _full_run(tmp_path)
    seq = make_plastome(**kw)
    (out / "reference.fasta").write_text(f">c1 a plastid\n{seq}\n")
    rows = "".join(f"c1\t{pos}\t.\t{seq[pos - 1]}\t{'A' if seq[pos - 1] != 'A' else 'C'}\t.\t.\t.\tGT\t1\t0\n"
                   for pos in positions)
    (out / "4_compared" / "ska" / "snps.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n" + rows)
    if features is not None:
        (out / "annotation.gb").write_text(_plastid_genbank(seq, features))
        info = json.loads((out / "run_info.json").read_text())
        info["annotation"] = {"file": "/x/NC_1.gb", "format": "genbank", "copy": "annotation.gb",
                              "transl_tables": [11]}
        (out / "run_info.json").write_text(json.dumps(info))
    return out, seq


DETECTED = "the inverted repeat detected in the reference sequence (two copies of 6,000 bp, 100% identical)"


def test_report_without_annotation_draws_the_regions_detected_in_the_sequence(tmp_path):
    out, _ = _plastid_run(tmp_path)
    page = build_report(out)
    assert 'class="band-ir"' in page and "LSC 9,000 bp" in page and "IRb 6,000 bp" in page
    assert "SSC 4,000 bp" in page and "IRa 6,000 bp" in page
    assert f"The band under the axis shows the LSC/IRb/SSC/IRa regions from {DETECTED}." in page
    assert "Hover a tick for the position, the alleles and the number of genomes with the alternate allele and " \
           "the region." in page
    assert "called in every genome. IRb</title>" in page and "called in every genome. LSC</title>" in page
    assert 'class="gene-' not in page and "<h3>SNPs</h3>" not in page and "Genes from" not in page  # No annotation
    assert f"The LSC/IRb/SSC/IRa regions were derived from {DETECTED}. The inverted repeat is found by matching " \
           "the k-mers" in page  # Methods
    # A reference without a repeat: no band, nothing said
    (out / "reference.fasta").write_text(">c1 a plastid\n" + "ACGT" * 6250 + "\n")
    page = build_report(out)
    assert 'class="band-ir"' not in page and "LSC/IRb/SSC/IRa" not in page


def test_report_with_an_annotation_without_repeat_features_detects_them(tmp_path):
    out, _ = _plastid_run(tmp_path, features=(
        "     gene            101..400\n                     /gene=\"psbA\"\n"
        "     CDS             101..400\n                     /gene=\"psbA\"\n"))
    page = build_report(out)
    assert 'class="band-ir"' in page and 'class="gene-cds"' in page
    assert f"The band above the genes shows the LSC/IRb/SSC/IRa regions from {DETECTED}." in page
    assert "2 SNPs on the annotated sequences: 1 in the LSC, 1 in IRb" in page  # The SNP table has the regions
    assert '<td class="" data-v="IRb">IRb</td>' in page
    assert "Genes were read from the annotation NC_1.gb" in page
    assert f"The LSC/IRb/SSC/IRa regions were derived from {DETECTED}." in page
    # With the repeats annotated (even a little off), the annotation's coordinates are kept and said to be used
    (out / "annotation.gb").write_text(_plastid_genbank(_plastid_run(tmp_path / "again")[1], (
        "     repeat_region   9001..15000\n                     /rpt_type=inverted\n"
        "     repeat_region   19001..24990\n                     /rpt_type=inverted\n")))
    page = build_report(out)
    assert "IRa 5,990 bp" in page and "regions from the annotated inverted repeats." in page
    assert "The LSC/IRb/SSC/IRa regions were derived from the annotated inverted repeats." in page
    assert "detected" not in page


def test_report_labels_genes_without_a_symbol_by_locus_tag_and_product(tmp_path):
    out, _ = _plastid_run(tmp_path, positions=(150, 250, 500, 12000), features=(
        "     gene            101..400\n                     /locus_tag=\"LK299_pgp087\"\n"
        "     CDS             101..400\n                     /locus_tag=\"LK299_pgp087\"\n"
        "                     /product=\"maturase K\"\n"
        "     gene            601..900\n                     /locus_tag=\"SIM_p010\"\n"
        "     CDS             601..900\n                     /locus_tag=\"SIM_p010\"\n"
        "                     /product=\"ribosomal <protein> & S12\"\n"))
    page = build_report(out)
    assert '<td class="gene" data-v="LK299_pgp087 (maturase K)">LK299_pgp087 (maturase K)</td>' in page  # Table
    assert "Genes with the most SNPs: <i>LK299_pgp087</i> (maturase K; 2)." in page  # Summary
    assert 'class="t-tiny t-ink2 t-gene">maturase K<title>LK299_pgp087 (maturase K): protein-coding gene' in page
    assert "intergenic between LK299_pgp087 (maturase K) and SIM_p010 (ribosomal &lt;protein&gt; &amp; S12)" \
        in page  # The context of the SNP at 500, escaped
    assert "ribosomal <protein>" not in page
    assert "LK299_pgp087 (maturase K) (CDS)" in page  # The tick's hover
    assert "SIM_p010 (ribosomal &lt;protein&gt; &amp; S12): protein-coding gene, + strand, 601–900 (300 bp)" \
        in page and "(300 bp), ribosomal" not in page  # The product is not repeated after the label


def test_methods_text_describes_the_regions_only_when_the_report_shows_a_band():
    from bacon.annotation import Region, RegionBand
    from bacon.report import methods_text
    base = {"settings": {}, "reference": {"file": "/x/ref.fa"}}
    assert "LSC/IRb/SSC/IRa" not in methods_text(base)
    # The regions recorded in run_info.json are not described when no map (no band) is in the report
    detected = {"source": "sequence", "regions": [{"name": "LSC", "start": 1, "end": 9000, "length": 9000}],
                "repeat": {"copies": [[9001, 15000], [19001, 25000]], "lengths": [6000, 6000], "differences": 0,
                           "identity": 1.0}}
    info = {**base, "reference": {"file": "/x/ref.fa", "regions": {"c1": detected, "c2": "none"}}}
    assert "LSC/IRb/SSC/IRa" not in methods_text(info) and "LSC/IRb/SSC/IRa" not in methods_text(info, None, {})
    # The bands of the map are described, one text when they all have the same source, else per sequence
    band = RegionBand([Region("LSC", 1, 9000, 9000)], "annotation")
    assert "derived from the annotated inverted repeats." in methods_text(info, None, {"c1": band})
    from bacon.annotation import InvertedRepeat
    found = RegionBand([Region("LSC", 1, 9000, 9000)], "sequence",
                       InvertedRepeat((9001, 15000), (19001, 25000), (6000, 6000), 0))
    text = methods_text(info, None, {"c1": found, "c2": band, "c3": RegionBand([], "none")})
    assert f"derived from c1: {DETECTED}; c2: the annotated inverted repeats." in text
    assert "k-mers" in text and "c3" not in text  # How the repeat is found, when one was; no band: not mentioned


def test_dense_gene_titles_use_the_labels(tmp_path, monkeypatch):
    import bacon.report as report
    from bacon.annotation import load_annotation
    gb = tmp_path / "a.gb"
    gb.write_text(_plastid_genbank("ACGT" * 1000, (
        "     gene            1..300\n                     /locus_tag=\"LK299_pgr007\"\n"
        "     rRNA            1..300\n                     /locus_tag=\"LK299_pgr007\"\n"
        "                     /product=\"23S ribosomal RNA\"\n"
        "     gene            complement(2000..2072)\n                     /locus_tag=\"OrsajCt141\"\n"
        "     tRNA            complement(2000..2072)\n                     /locus_tag=\"OrsajCt141\"\n"
        "                     /product=\"tRNA-Val\"\n")))
    annotation = load_annotation(gb, [("c1", 4000)])
    monkeypatch.setattr(report, "MAX_GENE_RECTS", 1)
    out = genome_map([("c1", 4000)], [], 2, annotation=annotation)
    assert "1 gene: LK299_pgr007 (23S ribosomal RNA)" in out and "1 gene: OrsajCt141 (tRNA-Val)" in out


def test_label_widths_follow_the_map_label():
    from bacon.annotation import Gene
    from bacon.report import _text_px, place_labels
    left, right = 118, 1100
    scale = (right - left) / 4000

    def gene(name: str, product: str, start: int) -> Gene:
        g = Gene("c1", name, "rRNA", 1, [(start, start + 100)], product=product)
        g.snps = 2
        return g

    short, long_ = _text_px("tRNA-Val", 9.5), _text_px("LK299_pgr007", 9.5)
    assert short + 1 < long_
    apart = int((short + 7) / scale)  # Two labels this far apart fit side by side when short, not when long
    placed, rows = place_labels([gene("LK299_pgr007", "tRNA-Val", 1000), gene("LK299_pgr008", "tRNA-Ile", 1000 + apart)],
                                left, right, scale)
    assert rows == 1 and [row for _, row, _ in placed] == [0, 0]
    placed, rows = place_labels([gene("LK299_pgr007", "23S ribosomal RNA", 1000),
                                 gene("LK299_pgr008", "16S ribosomal RNA", 1000 + apart)], left, right, scale)
    assert rows == 2 and [row for _, row, _ in placed] == [0, 1]  # Labelled by their 12-character identifiers


def test_band_height_is_reserved_without_annotation():
    import re

    from bacon.annotation import find_regions
    from tests.test_annotation import make_plastome
    seq = make_plastome()
    band = find_regions([], len(seq), seq)
    plain = genome_map([("c1", len(seq))], [], 2)
    banded = genome_map([("c1", len(seq))], [], 2, bands={"c1": band})

    def snp_row(svg: str) -> float:
        return float(re.search(r'y="([\d.]+)"[^>]*>SNPs</text>', svg).group(1))

    def height(svg: str) -> float:
        return float(re.search(r'viewBox="0 0 \d+ (\d+)"', svg).group(1))

    assert 'class="band-ir"' in banded and 'class="band-ir"' not in plain
    assert snp_row(banded) - snp_row(plain) == 16  # The SNP track moves down by the band's height ...
    assert height(banded) - height(plain) == 16  # ... and the figure is taller by as much


def test_merged_ticks_hover_gives_the_region_without_annotation():
    from bacon.annotation import find_regions
    from bacon.report import Snp
    from tests.test_annotation import make_plastome
    seq = make_plastome()
    band = find_regions([], len(seq), seq)
    snps = [Snp("c1", 100, "A", "G", 1, 0), Snp("c1", 101, "C", "T", 2, 0), Snp("c1", 12000, "A", "G", 1, 0)]
    out = genome_map([("c1", len(seq))], snps, 3, bands={"c1": band})
    assert "2 SNPs at c1:100–101; 0 with a missing call. 100 A&gt;G (1 alt) LSC 101 C&gt;T (2 alt) LSC</title>" in out
    assert "called in every genome. IRb</title>" in out


def test_caption_says_when_a_sequence_is_too_long_to_search(tmp_path):
    out = _full_run(tmp_path)
    (out / "reference.fasta").write_text(">c1 big\n" + "A" * 2_000_001 + "\n")
    page = build_report(out)
    assert "No search for an inverted repeat in sequences longer than 2 Mb." in page
    assert "LSC/IRb/SSC/IRa" not in page


def test_caption_places_the_band_per_sequence(tmp_path):
    from tests.test_annotation import make_plastome
    out = _full_run(tmp_path)
    s1, s2 = make_plastome(seed=1), make_plastome(seed=2)
    (out / "reference.fasta").write_text(f">c1 a\n{s1}\n>c2 b\n{s2}\n")
    (out / "4_compared" / "ska" / "snps.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n"
        f"c1\t150\t.\t{s1[149]}\t{'A' if s1[149] != 'A' else 'C'}\t.\t.\t.\tGT\t1\t0\n"
        f"c2\t12000\t.\t{s2[11999]}\t{'A' if s2[11999] != 'A' else 'C'}\t.\t.\t.\tGT\t1\t0\n")
    (out / "annotation.gb").write_text(_plastid_genbank(s1, (  # c1 only
        "     gene            101..400\n                     /gene=\"psbA\"\n"
        "     CDS             101..400\n                     /gene=\"psbA\"\n")))
    info = json.loads((out / "run_info.json").read_text())
    info["annotation"] = {"file": "/x/NC_1.gb", "format": "genbank", "copy": "annotation.gb", "transl_tables": [11]}
    (out / "run_info.json").write_text(json.dumps(info))
    page = build_report(out)
    assert f"The band above the genes, or under the axis for sequences without annotation shows the " \
           f"LSC/IRb/SSC/IRa regions from {DETECTED}." in page
    assert "with the gene, its context and the effect of the SNP, or the region alone on sequences without " \
           "annotation." in page
    assert "called in every genome. IRb</title>" in page  # c2's tick: the region alone
