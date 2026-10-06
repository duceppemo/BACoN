import json

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
    html_text = _genome_map_section(out, vcf, reference, {"assembler": "samtools"}, _Figures(), 1)
    assert "none of the 1 templated assemblies has an N base" in html_text and "N per" not in html_text
