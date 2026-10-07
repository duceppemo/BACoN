"""Reference annotations: GenBank and GFF3 parsing, regions, and the effects of SNPs on coding sequences."""

import gzip
import textwrap

import pytest

from bacon import BaconError
from bacon.annotation import (
    Cds,
    annotate_snps,
    annotation_format,
    cds_effect,
    derive_regions,
    genbank_fasta_records,
    genes_from_genbank,
    genes_from_gff3,
    load_annotation,
    parse_location,
    read_genbank,
    read_gff3,
    translate,
)

# A 120 bp sequence with, on the + strand, a split CDS (abc: two exons around an intron) and, on the - strand, a
# tRNA (trnX) and a plain CDS (rev); the rest is spacer.
PIECES = [("ATGGCTAAA", "abc exon 1 (1..9): M A K"), ("GTTAGCCCTTTA", "intron (10..21)"),
          ("GGTCGATAA", "abc exon 2 (22..30): G R *"), ("ACCCGGGTTT", "spacer (31..40)"),
          ("AAACCCGGGTTTAAATTTGG", "trnX, complement (41..60)"), ("GAAACCCAT", "spacer (61..69)"),
          ("CTATTTAAAGGCCAT", "rev, complement (70..84): ATG GCC TTT AAA TAG = M A F K *"),
          ("GGCCCTTAAATTTCCCGGGAAATTTGGGCCCAAATT", "spacer (85..120)")]
SEQ = "".join(piece for piece, _ in PIECES)
assert len(SEQ) == 120


def _gb(features: str, seq: str = SEQ, name: str = "ref", version: str = "ref.1", definition: str = "Test record.",
        circular: bool = True) -> str:
    origin = "\n".join(f"{i + 1:>9} " + " ".join(seq[j:j + 10].lower() for j in range(i, min(i + 60, len(seq)), 10))
                       for i in range(0, len(seq), 60))
    return (f"LOCUS       {name}  {len(seq)} bp    DNA     {'circular' if circular else 'linear'} PLN 01-JAN-2026\n"
            f"DEFINITION  {definition}\n"
            f"ACCESSION   {name}\n"
            + (f"VERSION     {version}\n" if version else "")
            + "FEATURES             Location/Qualifiers\n"
            + textwrap.indent(textwrap.dedent(features).strip("\n"), "     ") + "\n"
            + "ORIGIN      \n" + origin + "\n//\n")


FEATURES = """
gene            1..32
                /gene="abc"
                /locus_tag="T_001"
CDS             join(1..9,22..32)
                /gene="abc"
                /locus_tag="T_001"
                /codon_start=1
                /transl_table=11
                /product="protein ABC, a long
                product name"
                /translation="MAKGR"
gene            complement(41..60)
                /gene="trnX"
tRNA            complement(41..60)
                /gene="trnX"
                /product="tRNA-Xaa"
gene            complement(70..84)
                /gene="rev"
CDS             complement(70..84)
                /gene="rev"
                /codon_start=1
                /transl_table=11
                /product="reverse protein"
misc_feature    61..69
                /note="spacer"
"""

GFF = """\
##gff-version 3
##sequence-region ref.1 1 120
ref.1\ttest\tregion\t1\t120\t.\t+\t.\tID=ref.1;Is_circular=true
ref.1\ttest\tgene\t1\t32\t.\t+\t.\tID=gene-abc;Name=abc;gene=abc;locus_tag=T_001
ref.1\ttest\tmRNA\t1\t32\t.\t+\t.\tID=rna-abc;Parent=gene-abc
ref.1\ttest\tCDS\t1\t9\t.\t+\t0\tID=cds-abc;Parent=rna-abc;gene=abc;product=protein ABC%2C a long product name;transl_table=11
ref.1\ttest\tCDS\t22\t32\t.\t+\t0\tID=cds-abc;Parent=rna-abc;gene=abc;product=protein ABC%2C a long product name;transl_table=11
ref.1\ttest\tgene\t41\t60\t.\t-\t.\tID=gene-trnX;Name=trnX;gene=trnX
ref.1\ttest\ttRNA\t41\t60\t.\t-\t.\tID=rna-trnX;Parent=gene-trnX;product=tRNA-Xaa
ref.1\ttest\texon\t41\t60\t.\t-\t.\tID=exon-trnX;Parent=rna-trnX
ref.1\ttest\tgene\t70\t84\t.\t-\t.\tID=gene-rev;Name=rev;gene=rev
ref.1\ttest\tCDS\t70\t84\t.\t-\t0\tID=cds-rev;Parent=gene-rev;gene=rev;product=reverse protein
"""


# ---------------------------------------------------------------------------------------------------------------
# Locations
# ---------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("text, strand, parts, partial, ordered", [
    ("123..456", 1, [(123, 456)], (False, False), None),
    ("complement(123..456)", -1, [(123, 456)], (False, False), None),
    ("join(1..9,22..32)", 1, [(1, 9), (22, 32)], (False, False), None),
    ("order(1..9,22..32)", 1, [(1, 9), (22, 32)], (False, False), None),
    ("complement(join(10..20,30..40))", -1, [(10, 20), (30, 40)], (False, False), [(30, 40), (10, 20)]),
    ("join(complement(30..40),complement(10..20))", -1, [(10, 20), (30, 40)], (False, False), [(30, 40), (10, 20)]),
    ("complement(join(99091..99116,99653..99884,71485..71598))", -1,
     [(71485, 71598), (99091, 99116), (99653, 99884)], (False, False),
     [(71485, 71598), (99653, 99884), (99091, 99116)]),  # Trans-spliced rps12: exon 1 is the last listed
    ("<1..>100", 1, [(1, 100)], (True, True), None),
    ("complement(<5..80)", -1, [(5, 80)], (True, False), None),
    ("join(<1..9,22..>32)", 1, [(1, 9), (22, 32)], (True, True), None),
    ("42", 1, [(42, 42)], (False, False), None),
    ("12^13", 1, [(12, 13)], (False, False), None),
    ("join(J00194.1:100..202,1..10)", 1, [(1, 10)], (False, False), None),  # The remote part is skipped
    ("join( 1..9 , 22..32 )", 1, [(1, 9), (22, 32)], (False, False), None),
])
def test_parse_location(text, strand, parts, partial, ordered):
    loc = parse_location(text)
    assert (loc.strand, loc.parts, (loc.partial_low, loc.partial_high)) == (strand, parts, partial)
    assert [(s, e) for s, e, _ in loc.ordered] == (ordered or (parts[::-1] if strand < 0 else parts))
    assert all(st == strand for _, _, st in loc.ordered)


def test_parse_location_remote_only_and_errors():
    assert parse_location("J00194.1:100..202") is None
    for bad in ("", "abc", "10..5", "0..5", "join(1..5", "complement()"):
        with pytest.raises(ValueError):
            parse_location(bad)


# ---------------------------------------------------------------------------------------------------------------
# GenBank and GFF3 to the same genes
# ---------------------------------------------------------------------------------------------------------------

def _same_genes(genes):
    return [(g.name, g.kind, g.strand, g.extent, g.exons,
             [(c.parts, c.codon_start, c.transl_table) for c in g.cds], g.product) for g in genes]


def test_genbank_and_gff3_give_the_same_genes(tmp_path):
    gb = tmp_path / "a.gb"
    gb.write_text(_gb(FEATURES))
    records = read_genbank(gb)
    assert len(records) == 1
    rec = records[0]
    assert (rec.name, rec.definition, rec.length, rec.circular) == ("ref.1", "Test record", 120, True)
    assert rec.header == "ref.1 Test record" and rec.seq.upper() == SEQ
    cds = [f for f in rec.features if f.type == "CDS"][0]
    assert cds.qualifiers["product"] == "protein ABC, a long product name"  # Multi-line qualifier
    assert cds.qualifiers["translation"] == "MAKGR" and cds.qualifiers["codon_start"] == "1"
    from_gb = genes_from_genbank(rec.features)
    gff = tmp_path / "a.gff3"
    gff.write_text(GFF)
    features, lengths, circular = read_gff3(gff)
    assert lengths == {"ref.1": 120} and circular == {"ref.1"}
    from_gff = genes_from_gff3(features)
    assert _same_genes(from_gb) == _same_genes(from_gff) == [
        ("abc", "CDS", 1, [(1, 32)], [(1, 9), (22, 32)], [([(1, 9), (22, 32)], 1, 11)],
         "protein ABC, a long product name"),
        ("trnX", "tRNA", -1, [(41, 60)], [(41, 60)], [], "tRNA-Xaa"),
        ("rev", "CDS", -1, [(70, 84)], [(70, 84)], [([(70, 84)], 1, 11)], "reverse protein"),
    ]


def test_genbank_cds_without_gene_feature_pseudo_and_duplicate_names():
    rec = read_genbank_text(_gb("""
    CDS             10..30
                    /locus_tag="L1"
                    /pseudo
    gene            40..50
                    /gene="dup"
    tRNA            40..50
                    /gene="dup"
    gene            60..70
                    /gene="dup"
    rRNA            60..70
                    /gene="dup"
    gene            80..90
                    /gene="lonely"
    """))
    genes = genes_from_genbank(rec.features)
    assert [(g.name, g.kind, g.extent) for g in genes] == [
        ("L1", "pseudogene", [(10, 30)]), ("dup", "tRNA", [(40, 50)]), ("dup", "rRNA", [(60, 70)]),
        ("lonely", "other", [(80, 90)])]


def read_genbank_text(text: str):
    import io

    from bacon.annotation import _genbank_records
    return _genbank_records(io.StringIO(text))[0]


def test_gff3_pseudogene_phase_and_parentless_cds(tmp_path):
    gff = tmp_path / "p.gff3"
    gff.write_text("##gff-version 3\n"
                   "c\tx\tpseudogene\t1\t30\t.\t+\t.\tID=p1;gene=psA\n"
                   "c\tx\tCDS\t1\t30\t.\t+\t0\tID=cp1;Parent=p1;pseudo=true\n"
                   "c\tx\tCDS\t40\t60\t.\t-\t2\tID=c2;product=orphan\n"
                   "c\tx\tgene\t70\t90\t.\t+\t.\tID=g3;gene_biotype=pseudogene;Name=g3\n"
                   "c\tx\tgene\t100\t120\t.\t+\t.\tID=gene:E1;biotype=processed_pseudogene;Name=e1\n"  # Ensembl
                   "c\tx\tgene\t130\t160\t.\t+\t.\tID=g5;Name=noid\n"
                   "c\tx\tCDS\t130\t141\t.\t+\t0\tParent=g5\n"  # CDS lines without an ID: one CDS
                   "c\tx\tCDS\t150\t160\t.\t+\t0\tParent=g5\n")
    features, _, _ = read_gff3(gff)
    genes = genes_from_gff3(features)
    assert [(g.name, g.kind, g.strand) for g in genes] == [
        ("psA", "pseudogene", 1), ("orphan", "CDS", -1), ("g3", "pseudogene", 1), ("e1", "pseudogene", 1),
        ("noid", "CDS", 1)]
    assert genes[1].cds[0].codon_start == 3  # Phase 2 on the first part
    assert len(genes[4].cds) == 1 and genes[4].cds[0].parts == [(130, 141), (150, 160)]


def test_gff3_parts_in_the_order_of_translation_like_genbank(tmp_path):
    # The same coding sequences as a GenBank location and as GFF3 lines (start, end, strand, phase[, part], in
    # file order): a CDS trans-spliced across strands (a), three reverse-strand parts listed 5' to 3' out of
    # coordinate order like the IRb copy of rps12 (b), a CDS across the origin (c), a reverse-strand spliced CDS
    # listed by ascending coordinate (d, Ensembl) and 5' to 3' (e, NCBI), and one across the origin on the
    # reverse strand with NCBI's `part=` numbers (f) or listed by ascending coordinate without them (g)
    cases = {
        "a": ("join(complement(10..18),30..38)", [(10, 18, "-", 0), (30, 38, "+", 0)]),
        "b": ("complement(join(50..55,70..78,40..45))", [(40, 45, "-", 0), (70, 78, "-", 0), (50, 55, "-", 0)]),
        "c": ("join(110..120,1..9)", [(110, 120, "+", 0), (1, 9, "+", 1)]),
        "d": ("complement(join(70..75,77..84))", [(70, 75, "-", 1), (77, 84, "-", 0)]),
        "e": ("complement(join(70..75,77..84))", [(77, 84, "-", 0), (70, 75, "-", 1)]),
        "f": ("complement(join(110..120,1..9))", [(110, 120, "-", 0, 2), (1, 9, "-", 0, 1)]),
        "g": ("complement(join(110..120,1..9))", [(1, 9, "-", 0), (110, 120, "-", 0)]),  # f, Ensembl order
    }
    gb = tmp_path / "a.gb"
    gb.write_text(_gb("".join(f'CDS             {loc}\n                /gene="{name}"\n'
                              for name, (loc, _) in cases.items())))
    gff = tmp_path / "a.gff3"
    gff.write_text("##gff-version 3\n##sequence-region ref.1 1 120\n" + "".join(
        f"ref.1\tt\tCDS\t{s}\t{e}\t.\t{strand}\t{phase}\tID=cds-{name};gene={name}"
        + (f";part={rest[0]}" if rest else "") + "\n"
        for name, (_, lines) in cases.items() for s, e, strand, phase, *rest in lines))
    anns = [load_annotation(path, [("ref.1", 120)]) for path in (gb, gff)]
    from_gb, from_gff = ({g.name: g for g in ann.sequences["ref.1"].genes} for ann in anns)

    def same(g):
        return g.extent, g.strand, [(c.parts, c.strands, c.codon_start, c.coding_sequence(SEQ)) for c in g.cds]

    for name in cases:
        assert same(from_gff[name]) == same(from_gb[name]), name
    assert from_gb["a"].cds[0].parts == [(10, 18), (30, 38)] and from_gb["a"].cds[0].strands == [-1, 1]
    assert from_gb["b"].cds[0].parts == [(40, 45), (70, 78), (50, 55)] and from_gb["b"].strand == -1
    assert from_gb["c"].cds[0].parts == [(110, 120), (1, 9)] and from_gb["c"].extent == [(1, 9), (110, 120)]
    assert from_gb["c"].cds[0].coding_sequence(SEQ) == SEQ[109:] + SEQ[:9]
    assert from_gb["d"].cds[0].parts == [(77, 84), (70, 75)] == from_gb["e"].cds[0].parts
    assert from_gb["f"].cds[0].parts == [(1, 9), (110, 120)] and from_gb["f"].strand == -1
    for pos in (1, 5, 12, 42, 52, 72, 80, 115):
        effects = [[(e.gene, e.codons, e.change, e.kind) for e in _effects(ann, pos, SEQ[pos - 1], "A")[0].effects]
                   for ann in anns]
        assert effects[0] == effects[1] and effects[0], pos  # The same effects from both files
    # Across the origin: GGG CCC AAA TTA TGG... on the + strand, TTT AGC CAT AAT TTG... on the - strand
    assert _effects(anns[0], 1, "A", "G")[1] == _effects(anns[1], 1, "A", "G")[1] == [
        ("c", "TTA>TTG", "L4L", "synonymous"), ("f", "CAT>CAC", "H3H", "synonymous"),
        ("g", "CAT>CAC", "H3H", "synonymous")]


def test_feature_across_the_origin_without_gene_feature(tmp_path):
    features = """
    CDS             join(110..120,1..9)
                    /gene="wrap"
    tRNA            complement(join(100..120,1..3))
                    /product="tRNA-Wrap"
    CDS             join(20..30,70..80)
                    /gene="wide"
    """
    ann = _ann(tmp_path, features)
    genes = {g.name: g for g in ann.sequences["ref.1"].genes}
    assert genes["wrap"].extent == [(1, 9), (110, 120)] and genes["wrap"].length == 20
    assert genes["tRNA-Wrap"].extent == [(1, 3), (100, 120)] and genes["tRNA-Wrap"].strand == -1
    assert genes["wide"].extent == [(20, 30), (70, 80)]  # Spanning more than half the sequence: across the origin
    info, _ = _effects(ann, 50, SEQ[49], "A")
    assert not info.genes and info.context.startswith("intergenic")
    assert _effects(ann, 5, "C", "A")[0].context == "CDS" and _effects(ann, 115, SEQ[114], "A")[0].context == "CDS / tRNA"
    # Without the sequence's length, only parts out of order tell: the wide CDS is then one stretch
    genes = {g.name: g for g in genes_from_genbank(read_genbank_text(_gb(features)).features)}
    assert genes["wide"].extent == [(20, 80)] and genes["wrap"].extent == [(1, 9), (110, 120)]


def test_effect_through_a_shared_exon_given_once(tmp_path):
    # Like the two rps12 genes of a plastome, whose 5' exon is one and the same stretch of the LSC
    features = """
    gene            join(1..9,22..32)
                    /gene="rps12"
                    /locus_tag="L1"
    CDS             join(1..9,22..32)
                    /gene="rps12"
                    /locus_tag="L1"
    gene            join(1..9,36..46)
                    /gene="rps12"
                    /locus_tag="L2"
    CDS             join(1..9,36..46)
                    /gene="rps12"
                    /locus_tag="L2"
    """
    ann = _ann(tmp_path, features)
    info, effects = _effects(ann, 5, "C", "A")
    assert [g.name for g in info.genes] == ["rps12", "rps12"] and info.context == "CDS"
    assert effects == [("rps12", "GCT>GAT", "A2D", "missense")]
    assert _effects(ann, 23, "G", "A")[1] == [("rps12", "GGT>GAT", "G4D", "missense")]  # Exon 2 of L1 only
    assert _effects(ann, 5, "C", "A,T")[1] == [("rps12", "GCT>GAT", "A2D", "missense"),
                                              ("rps12", "GCT>GTT", "A2V", "missense")]


def test_transl_except_codon_gives_no_effect(tmp_path):
    from bacon.annotation import _transl_except
    seq = "ATGGCTTGAAAATAG" + "C" * 105  # M A Sec K *
    features = """
    CDS             1..15
                    /gene="sel"
                    /transl_except=(pos:7..9,aa:Sec)
    """
    ann = _ann(tmp_path, features, seq)
    assert ann.sequences["ref.1"].genes[0].cds[0].transl_except == [(7, 9)]
    for pos, ref, alt in ((7, "T", "C"), (8, "G", "A"), (9, "A", "G")):
        info, effects = _effects(ann, pos, ref, alt, seq)
        assert info.context == "CDS" and effects == [], pos  # Neither stop lost nor stop retained
    assert _effects(ann, 5, "C", "A", seq)[1] == [("sel", "GCT>GAT", "A2D", "missense")]
    assert _effects(ann, 11, "A", "G", seq)[1] == [("sel", "AAA>AGA", "K4R", "missense")]
    assert _transl_except("(pos:7..9,aa:Sec); (pos:complement(20..22),aa:TERM)") == [(7, 9), (20, 22)]
    assert _transl_except("(pos:4299366,aa:TERM)") == [(4299366, 4299366)] and _transl_except("") == []
    gff = tmp_path / "s.gff3"
    gff.write_text("##gff-version 3\n##sequence-region ref.1 1 120\n"
                   "ref.1\tt\tCDS\t1\t15\t.\t+\t0\tID=c;gene=sel;transl_except=(pos:7..9%2Caa:Sec)\n")
    ann = load_annotation(gff, [("ref.1", 120)])
    assert ann.sequences["ref.1"].genes[0].cds[0].transl_except == [(7, 9)]
    assert _effects(ann, 8, "G", "A", seq)[1] == []


def test_genbank_doubled_quotes():
    rec = read_genbank_text(_gb('''
    CDS             1..9
                    /gene="q"
                    /note="a ""quoted"" word; the line ends with ""
                    still the note"
                    /product="say ""hi"""
                    /standard_name=""
                    /db_xref="GeneID:1"
    '''))
    q = rec.features[0].qualifiers
    assert q["note"] == 'a "quoted" word; the line ends with " still the note'
    assert q["product"] == 'say "hi"' and q["standard_name"] == "" and q["db_xref"] == "GeneID:1"


def test_gzipped_inputs_and_format_detection(tmp_path):
    gb = tmp_path / "a.gbk.gz"
    with gzip.open(gb, "wt") as fh:
        fh.write(_gb(FEATURES))
    assert annotation_format(gb) == "genbank"
    assert len(read_genbank(gb)[0].features) == 7
    gff = tmp_path / "a.gff.gz"
    with gzip.open(gff, "wt") as fh:
        fh.write(GFF)
    assert annotation_format(gff) == "gff3"
    assert len(read_gff3(gff)[0]) == 9
    by_content = tmp_path / "noext"
    by_content.write_text(_gb(FEATURES))
    assert annotation_format(by_content) == "genbank"
    by_content.write_text(GFF)
    assert annotation_format(by_content) == "gff3"
    by_content.write_text(">fasta\nACGT\n")
    assert annotation_format(by_content) is None
    assert [r.header for r in genbank_fasta_records(gb)] == ["ref.1 Test record"]
    assert genbank_fasta_records(gb)[0].seq.upper() == SEQ


def test_genbank_without_sequence_cannot_be_the_reference(tmp_path):
    gb = tmp_path / "a.gb"
    gb.write_text(_gb(FEATURES).split("ORIGIN")[0] + "//\n")
    with pytest.raises(BaconError, match="no sequence .*--annotation"):
        genbank_fasta_records(gb)
    assert load_annotation(gb, [("ref.1", 120)]).genes == 3  # Still a valid annotation


# ---------------------------------------------------------------------------------------------------------------
# Loading against a reference: names, lengths, warnings
# ---------------------------------------------------------------------------------------------------------------

def test_load_annotation_matches_names_and_warns(tmp_path):
    gb = tmp_path / "a.gb"
    gb.write_text(_gb(FEATURES))
    ann = load_annotation(gb, [("ref.1", 120), ("other", 50)])
    assert ann.format == "genbank" and ann.genes == 3 and list(ann.sequences) == ["ref.1"] and not ann.warnings
    assert ann.tables == [11] and not ann.has_regions
    # A single annotated sequence of the same length is taken for a single reference sequence of another name
    ann = load_annotation(gb, [("chr", 120)])
    assert list(ann.sequences) == ["chr"] and "taken for the reference sequence 'chr'" in ann.warnings[0]
    assert ann.sequences["chr"].genes[0].seq == "chr"
    # Different length: no match, a warning, no features
    ann = load_annotation(gb, [("chr", 121)])
    assert ann.sequences == {} and "match no reference sequence by name" in ann.warnings[0]
    # A different length under the same name is said; features beyond the end of the sequence are dropped
    ann = load_annotation(gb, [("ref.1", 50)])
    assert ann.genes == 1 and ann.skipped == 5 and len(ann.warnings) == 2
    assert "the annotation of ref.1 is 120 bp, the reference 50 bp: is it the annotation" in ann.warnings[0]
    assert "5 feature(s) beyond the end" in ann.warnings[1]
    gff = tmp_path / "a.gff3"
    gff.write_text(GFF)
    ann = load_annotation(gff, [("ref.1", 121)])  # The ##sequence-region length
    assert ann.genes == 3 and "the annotation of ref.1 is 120 bp, the reference 121 bp" in ann.warnings[0]


def test_load_annotation_errors(tmp_path):
    with pytest.raises(BaconError, match="not found"):
        load_annotation(tmp_path / "missing.gb", [("a", 1)])
    fasta = tmp_path / "x.fasta"
    fasta.write_text(">a\nACGT\n")
    with pytest.raises(BaconError, match="not a GenBank .* or GFF3"):
        load_annotation(fasta, [("a", 4)])
    empty = tmp_path / "e.gff3"
    empty.write_text("##gff-version 3\n")
    with pytest.raises(BaconError, match="no feature"):
        load_annotation(empty, [("a", 4)])
    bad = tmp_path / "bad.gb.gz"
    bad.write_bytes(b"\x1f\x8b\x08\x00garbage")
    with pytest.raises(BaconError, match="truncated or corrupt"):
        load_annotation(bad, [("a", 4)])


# ---------------------------------------------------------------------------------------------------------------
# Regions
# ---------------------------------------------------------------------------------------------------------------

def _features(text):
    return read_genbank_text(_gb(text, seq="A" * 10000)).features


def test_regions_from_named_inverted_repeats():
    features = _features("""
    misc_feature    1..4000
                    /note="LSC"
    misc_feature    4001..6000
                    /note="IRB"
    misc_feature    6001..7000
                    /note="SSC"
    misc_feature    7001..9000
                    /note="IRA"
    misc_feature    4000..4001
                    /note="JLB; junction LSC-IRB"
    """)
    regions = derive_regions(features, 10000)
    assert [(r.name, r.start, r.end, r.length) for r in regions] == [
        ("IRb", 4001, 6000, 2000), ("SSC", 6001, 7000, 1000), ("IRa", 7001, 9000, 2000), ("LSC", 9001, 4000, 5000)]
    lsc = regions[-1]
    assert lsc.contains(9500) and lsc.contains(5) and not lsc.contains(5000)


def test_regions_named_by_convention_and_wrap():
    # Unnamed inverted repeats: the larger gap is the LSC and the repeat after it is IRb; here the SSC lies
    # between the repeats and the LSC spans the origin.
    features = _features("""
    repeat_region   3001..5000
                    /rpt_type=inverted
    repeat_region   6001..8000
                    /rpt_type=inverted
    """)
    assert [(r.name, r.start, r.end, r.length) for r in derive_regions(features, 10000)] == [
        ("IRb", 3001, 5000, 2000), ("SSC", 5001, 6000, 1000), ("IRa", 6001, 8000, 2000), ("LSC", 8001, 3000, 5000)]
    # The sequence starts in the SSC: the gap between the repeats is the LSC
    features = _features("""
    repeat_region   1001..2000
                    /note="inverted repeat A"
    repeat_region   8001..9000
                    /note="inverted repeat B"
    """)
    assert [(r.name, r.start, r.end) for r in derive_regions(features, 10000)] == [
        ("IRa", 1001, 2000), ("LSC", 2001, 8000), ("IRb", 8001, 9000), ("SSC", 9001, 1000)]
    # A sequence starting exactly at a junction: no wrap
    features = _features("""
    repeat_region   1..2000
                    /rpt_type=inverted
    repeat_region   7001..9000
                    /rpt_type=inverted
    """)
    assert [(r.name, r.start, r.end) for r in derive_regions(features, 10000)] == [
        ("IRa", 1, 2000), ("LSC", 2001, 7000), ("IRb", 7001, 9000), ("SSC", 9001, 10000)]


def test_no_regions_without_two_inverted_repeats():
    assert derive_regions(_features("""
    repeat_region   10..40
                    /rpt_type=inverted
    repeat_region   900..930
                    /rpt_type=inverted
    """), 10000) == []  # Too short to be plastome repeats (hairpins)
    assert derive_regions(_features("""
    misc_feature    1..6000
                    /note="LSC"
    """), 10000) == []
    assert derive_regions(_features("""
    repeat_region   1..5000
                    /rpt_type=inverted
    repeat_region   5001..10000
                    /rpt_type=inverted
    """), 10000) == []  # No single-copy region between them


# ---------------------------------------------------------------------------------------------------------------
# Effects
# ---------------------------------------------------------------------------------------------------------------

def test_translation_tables():
    assert translate("ATG") == "M" and translate("TAA") == "*" and translate("TGG") == "W"
    assert translate("GTG", 11, start=True) == "M" and translate("GTG", 11) == "V"
    assert translate("GTG", 1, start=True) == "V" and translate("TGA", 4) == "W" and translate("TGA", 1) == "*"
    assert translate("ANN") == "X" and translate("TTA", 4, start=True) == "M"
    assert translate("GTG", 99, start=True) == "M"  # An unknown table: the standard code with table 11's starts
    # Mitochondrial codes (NCBI): distinctive codons and start codons
    assert [translate(c, 2) for c in ("AGA", "AGG", "ATA", "TGA", "AAA")] == ["*", "*", "M", "W", "K"]
    assert translate("ATT", 2, start=True) == "M" and translate("ATT", 2) == "I" and translate("TTG", 2, start=True) == "L"
    assert [translate(c, 3) for c in ("CTT", "CTC", "CTA", "CTG", "ATA", "TGA")] == ["T", "T", "T", "T", "M", "W"]
    assert translate("CTG", 3, start=True) == "T" and translate("GTG", 3, start=True) == "M"
    assert [translate(c, 5) for c in ("AGA", "AGG", "ATA", "TGA", "AAA")] == ["S", "S", "M", "W", "K"]
    assert translate("TTG", 5, start=True) == "M" and translate("ATC", 5, start=True) == "M"
    assert [translate(c, 9) for c in ("AAA", "AGA", "AGG", "TGA", "ATA")] == ["N", "S", "S", "W", "I"]
    assert translate("GTG", 9, start=True) == "M" and translate("ATA", 9, start=True) == "I"
    assert [translate(c, 13) for c in ("AGA", "AGG", "ATA", "TGA")] == ["G", "G", "M", "W"]
    assert translate("TTG", 13, start=True) == "M" and translate("ATT", 13, start=True) == "I"
    assert [translate(c, 14) for c in ("AAA", "AGA", "AGG", "TAA", "TGA", "TAG")] == ["N", "S", "S", "Y", "W", "*"]
    assert translate("GTG", 14, start=True) == "V" and translate("ATG", 14, start=True) == "M"
    assert translate("AGA", 4) == "R" and translate("AGA", 11) == "R"  # Unchanged in the other tables


def _ann(tmp_path, features=FEATURES, seq=SEQ):
    tmp_path.mkdir(exist_ok=True)
    gb = tmp_path / "a.gb"
    gb.write_text(_gb(features, seq=seq))
    return load_annotation(gb, [("ref.1", len(seq))])


def _effects(ann, pos, ref, alt, seq=SEQ):
    info = annotate_snps(ann, [("ref.1", pos, ref, alt)], {"ref.1": seq})[("ref.1", pos)]
    return info, [(e.gene, e.codons, e.change, e.kind) for e in info.effects]


def test_effects_forward_split_cds_and_introns(tmp_path):
    # abc: join(1..9,22..30): ATG GCT AAA | GGT CGA TAA -> M A K G R *
    ann = _ann(tmp_path)
    info, effects = _effects(ann, 5, "C", "A")  # GCT -> GAT
    assert info.context == "CDS" and info.region == "" and effects == [("abc", "GCT>GAT", "A2D", "missense")]
    assert _effects(ann, 6, "T", "C")[1] == [("abc", "GCT>GCC", "A2A", "synonymous")]
    assert _effects(ann, 7, "A", "T")[1] == [("abc", "AAA>TAA", "K3*", "nonsense")]
    assert _effects(ann, 22, "G", "A")[1] == [("abc", "GGT>AGT", "G4S", "missense")]  # First base of exon 2
    assert _effects(ann, 30, "A", "G")[1] == [("abc", "TAA>TAG", "*6*", "stop retained")]
    assert _effects(ann, 29, "A", "G")[1] == [("abc", "TAA>TGA", "*6*", "stop retained")]
    assert _effects(ann, 28, "T", "C")[1] == [("abc", "TAA>CAA", "*6Q", "stop lost")]
    assert _effects(ann, 1, "A", "G")[1] == [("abc", "ATG>GTG", "M1M", "start retained")]  # GTG starts (table 11)
    assert _effects(ann, 2, "T", "C")[1] == [("abc", "ATG>ACG", "M1T", "start lost")]
    info, effects = _effects(ann, 15, "G", "A")
    assert info.context == "intron" and effects == [] and [g.name for g in info.genes] == ["abc"]
    # With the standard code, GTG is not a start codon
    table1 = _ann(tmp_path / "t1", FEATURES.replace("/transl_table=11", "/transl_table=1"))
    assert _effects(table1, 1, "A", "G")[1] == [("abc", "ATG>GTG", "M1V", "start lost")]


def test_effects_reverse_strand_trna_intergenic_and_wrap(tmp_path):
    # rev: complement(70..84) reads ATG GCC TTT AAA TAG = M A F K *; the VCF alleles are on the + strand
    ann = _ann(tmp_path)
    assert ann.sequences["ref.1"].genes[2].cds[0].coding_sequence(SEQ) == "ATGGCCTTTAAATAG"
    assert SEQ[83] == "T" and SEQ[82] == "A" and SEQ[69] == "C"
    assert _effects(ann, 84, "T", "A")[1] == [("rev", "ATG>TTG", "M1M", "start retained")]  # Complement: A -> T
    assert _effects(ann, 83, "A", "T")[1] == [("rev", "ATG>AAG", "M1K", "start lost")]
    assert _effects(ann, 77, "A", "G")[1] == [("rev", "TTT>TCT", "F3S", "missense")]
    assert _effects(ann, 70, "C", "G")[1] == [("rev", "TAG>TAC", "*5Y", "stop lost")]
    info, effects = _effects(ann, 50, SEQ[49], "A")
    assert info.context == "tRNA" and effects == [] and [g.name for g in info.genes] == ["trnX"]
    info, effects = _effects(ann, 65, SEQ[64], "A")
    assert info.context == "intergenic between trnX and rev" and not info.genes
    assert _effects(ann, 100, SEQ[99], "A")[0].context == "intergenic between rev and abc"  # Circular: wraps
    assert _effects(ann, 35, SEQ[34], "A")[0].context == "intergenic between abc and trnX"
    linear = tmp_path / "linear"
    linear.mkdir()
    gb = linear / "a.gb"
    gb.write_text(_gb(FEATURES, circular=False))
    info = annotate_snps(load_annotation(gb, [("ref.1", 120)]), [("ref.1", 100, SEQ[99], "A")],
                         {"ref.1": SEQ})[("ref.1", 100)]
    assert info.context == "intergenic after rev"


def test_effects_codon_start_partial_multiallelic_and_overlap(tmp_path):
    features = """
    CDS             <1..12
                    /gene="p2"
                    /codon_start=2
    CDS             <1..13
                    /gene="p3"
                    /codon_start=3
    CDS             22..>35
                    /gene="tail"
    CDS             5..10
                    /gene="overlap"
    """
    ann = _ann(tmp_path, features)
    # codon_start=2: TGG CTA AAG G -> W L K (position 2 is the first base; partial, so not a start codon)
    assert _effects(ann, 2, "T", "A")[1] == [("p2", "TGG>AGG", "W1R", "missense")]
    # codon_start=3: GGC TAA AGG T -> G * R
    assert ("p3", "GGC>AGC", "G1S", "missense") in _effects(ann, 3, "G", "A")[1]
    # Position 1 is before the first codon of both: no effect, context CDS
    info, effects = _effects(ann, 1, "A", "G")
    assert info.context == "CDS" and effects == []
    # A partial 3' end: the incomplete last codon (34..35) gives no effect; the codon before does
    assert _effects(ann, 35, SEQ[34], "A")[1] == [] and _effects(ann, 34, SEQ[33], "A")[1] == []
    assert _effects(ann, 33, "C", "A")[1] == [("tail", "ACC>ACA", "T4T", "synonymous")]
    # Several alternate alleles: one effect per allele, in order
    assert _effects(ann, 2, "T", "A,C")[1] == [("p2", "TGG>AGG", "W1R", "missense"),
                                              ("p2", "TGG>CGG", "W1R", "missense")]
    # Three overlapping CDS: all reported (CTA is not a start codon, so overlap's first codon is L)
    assert _effects(ann, 7, "A", "T")[1] == [("p2", "CTA>CTT", "L2L", "synonymous"),
                                            ("p3", "TAA>TTA", "*2L", "stop lost"),
                                            ("overlap", "CTA>CTT", "L1L", "synonymous")]
    # A mismatching REF base gives no effect
    assert _effects(ann, 2, "C", "A")[1] == []


def test_trans_spliced_cds_across_strands(tmp_path):
    # Like the IRa copy of rps12: exon 1 on the - strand, exons 2 and 3 on the + strand, listed 5' to 3'
    features = """
    CDS             join(complement(10..18),30..38)
                    /gene="ts"
                    /trans_splicing
    """
    ann = _ann(tmp_path, features)
    gene = ann.sequences["ref.1"].genes[0]
    cds = gene.cds[0]
    assert gene.extent == [(10, 18), (30, 38)] and gene.strand == 1  # Most of its bases are on the + strand
    assert cds.parts == [(10, 18), (30, 38)] and cds.strands == [-1, 1]
    assert cds.coding_sequence(SEQ) == "AGGGCTAACAACCCGGGT"  # AGG GCT AAC AAC CCG GGT = R A N N P G
    assert _effects(ann, 18, "T", "C")[1] == [("ts", "AGG>GGG", "R1G", "missense")]  # Complemented part
    assert _effects(ann, 31, "A", "G")[1] == [("ts", "AAC>AGC", "N4S", "missense")]  # Forward part
    info, _ = _effects(ann, 25, SEQ[24], "A")
    assert info.context.startswith("intergenic") and not info.genes  # No intron between trans-spliced parts


def test_cds_effect_directly_and_n_in_codon():
    cds = Cds(1, [(1, 6)])
    assert cds_effect(cds, "g", 4, "A", "AT", "ATGAAA") is None  # Not a SNP
    assert cds_effect(cds, "g", 7, "A", "T", "ATGAAA") is None  # Outside
    assert cds_effect(cds, "g", 4, "A", "T", "ATGANA") is None  # N in the codon
    e = cds_effect(cds, "g", 4, "A", "T", "ATGAAA")
    assert (e.codons, e.change, e.kind, e.position) == ("AAA>TAA", "K2*", "nonsense", 2)


def test_trans_spliced_gene_lists_its_5_prime_exon_first(tmp_path):
    features = """
    gene            complement(join(60..80,10..20))
                    /gene="ts"
                    /trans_splicing
    CDS             complement(join(60..80,10..20))
                    /gene="ts"
                    /trans_splicing
    """
    ann = _ann(tmp_path, features)
    gene = ann.sequences["ref.1"].genes[0]
    assert gene.extent == [(10, 20), (60, 80)] and gene.cds[0].parts == [(10, 20), (60, 80)]  # 5' exon first
    assert gene.cds[0].strands is None and gene.strand == -1
    info, _ = _effects(ann, 40, SEQ[39], "A")
    assert info.context.startswith("intergenic") and not info.genes


def test_gene_snp_counts(tmp_path):
    ann = _ann(tmp_path)
    annotate_snps(ann, [("ref.1", 5, "C", "A"), ("ref.1", 6, "T", "C"), ("ref.1", 15, "G", "A"),
                        ("ref.1", 15, "G", "T"), ("other", 1, "A", "C")], {"ref.1": SEQ})
    assert [(g.name, g.snps) for g in ann.sequences["ref.1"].genes] == [("abc", 3), ("trnX", 0), ("rev", 0)]


def test_inverted_repeat_copies_named_in_free_text():
    from bacon.annotation import _IR_COPY
    for text, copy in [("IRa", "a"), ("IRB", "b"), ("inverted repeat B", "b"), ("inverted repeat IRb", "b"),
                       ("inverted repeat region IRa", "a"), ("Inverted Repeat A region", "a"), ("IR_A", "a"),
                       (" IRa (25,341 bp)", "a")]:
        match = _IR_COPY.match(text)
        assert match and match.group(1).lower() == copy, text
    for text in ("IRAK1 binding site", "inverted repeat", "spirAl", "LSC", "junction LSC-IRB",
                 "trans splicing 5'-rps12 and 3'-rps12 in IRA", "JLB"):
        assert not _IR_COPY.match(text), text


def test_regions_ignore_mentions_duplicates_and_other_inverted_repeats():
    # A note mentioning a copy does not annotate it; both repeats annotated twice (repeat_region and
    # misc_feature); a third inverted repeat (a transposon's) beside the two named ones; three unnamed ones
    both = """
    repeat_region   4001..6000
                    /rpt_type=inverted
                    /note="inverted repeat B"
    misc_feature    4001..6000
                    /note="IRB"
    repeat_region   7001..9000
                    /rpt_type=inverted
                    /note="inverted repeat A"
    misc_feature    7001..9000
                    /note="IRA"
    misc_feature    1001..2000
                    /note="trans splicing 5'-rps12 and 3'-rps12 in IRA"
    misc_feature    4000..4001
                    /note="JLB; junction LSC-IRB"
    """
    transposon = """
    repeat_region   100..700
                    /rpt_type=inverted
                    /note="terminal inverted repeat of a transposon"
    """
    unnamed = """
    repeat_region   4001..6000
                    /rpt_type=inverted
    repeat_region   7001..9000
                    /rpt_type=inverted
    """
    expected = [("IRb", 4001, 6000, 2000), ("SSC", 6001, 7000, 1000), ("IRa", 7001, 9000, 2000),
                ("LSC", 9001, 4000, 5000)]
    for text in (both, both + transposon, unnamed + transposon):
        assert [(r.name, r.start, r.end, r.length) for r in derive_regions(_features(text), 10000)] == expected


def test_inverted_repeat_across_the_origin():
    expected = [("SSC", 1701, 3000, 1300), ("IRa", 3001, 5000, 2000), ("LSC", 5001, 9700, 4700),
                ("IRb", 9701, 1700, 2000)]
    joined = """
    repeat_region   3001..5000
                    /rpt_type=inverted
    repeat_region   join(9701..10000,1..1700)
                    /rpt_type=inverted
    """
    in_two = """
    repeat_region   3001..5000
                    /note="IRa"
    repeat_region   9701..10000
                    /note="IRb"
    repeat_region   1..1700
                    /note="IRb"
    """
    for text in (joined, in_two):
        regions = derive_regions(_features(text), 10000)
        assert [(r.name, r.start, r.end, r.length) for r in regions] == expected
        assert regions[-1].contains(9800) and regions[-1].contains(100) and not regions[-1].contains(2000)


def test_5_prime_partial_cds_starting_with_atg_has_no_start_codon():
    cds = Cds(1, [(1, 12)], 1, 11, True, False)  # <1..12: the first codon is not the start
    e = cds_effect(cds, "g", 2, "T", "C", "ATGGCTAAATGA")
    assert (e.codons, e.change, e.kind) == ("ATG>ACG", "M1T", "missense")


def test_region_mention_in_a_wide_feature_note_does_not_label_it():
    from bacon.annotation import region_label
    repeats = """
    repeat_region   4001..4600
                    /rpt_type=inverted
                    /note="IRb"
    repeat_region   5401..6000
                    /rpt_type=inverted
                    /note="IRa"
    """
    wide = """
    misc_feature    1..9000
                    /note="trans splicing 5'-rps12 and 3'-rps12 in IRA"
    """
    [feature] = _features(wide)
    assert region_label(feature) is None

    def band(text):
        return [(r.name, r.start, r.end) for r in derive_regions(_features(text), 10000)]

    assert band(repeats + wide) == band(repeats) and sorted(n for n, _, _ in band(repeats)) == ["IRa", "IRb", "LSC", "SSC"]


def test_genbank_doubled_quotes_on_continuation_lines():
    rec = read_genbank_text(_gb('''
    CDS             1..9
                    /note="first line
                    a ""quoted"" word on the second, and one at its end: ""
                    then ""one"" closing the value"""
                    /gene="q"
    '''))
    q = rec.features[0].qualifiers
    assert q["note"] == 'first line a "quoted" word on the second, and one at its end: " then "one" closing the value"'
    assert q["gene"] == "q"


WRAPPED_CDS = "ATG" + "GCTAAAGAATTTCCCGGTACTCATGCTAAAGAATTTCCCGGTACTCATGCTAAA" + "TAA"  # 60 bp, 20 codons
WRAPPED_SEQ = WRAPPED_CDS[30:] + "A" * 60 + WRAPPED_CDS[:30]  # 120 bp: the CDS is 91..120 then 1..30 (+ strand)


def _wrapped_coding(gff, sequences: list[tuple[str, int]], strand: str) -> str:
    """The coding sequence read from the first gene of a GFF3 annotation, on the + strand of WRAPPED_SEQ or on
    the - strand of its reverse complement (where the same CDS is at 1..30 read downwards, then 91..120)."""
    from bacon.annotation import reverse_complement
    ann = load_annotation(gff, sequences)
    [sequence] = ann.sequences.values()
    [gene] = sequence.genes
    assert gene.extent == [(1, 30), (91, 120)]
    return gene.cds[0].coding_sequence(WRAPPED_SEQ if strand == "+" else reverse_complement(WRAPPED_SEQ))


def test_gff3_cds_across_the_origin_listed_by_coordinate(tmp_path):
    # A sorted GFF3 lists the parts of a CDS across the origin by ascending coordinate: the sequence length (from
    # ##sequence-region, else from the reference) tells them from a spliced CDS, on both strands
    for strand in "+-":
        lines = [f"chr\t.\tCDS\t1\t30\t.\t{strand}\t0\tID=c1;gene=wrap\n",
                 f"chr\t.\tCDS\t91\t120\t.\t{strand}\t0\tID=c1;gene=wrap\n"]
        gff = tmp_path / f"region{strand}.gff3"
        gff.write_text("##gff-version 3\n##sequence-region chr 1 120\n" + "".join(lines))
        assert _wrapped_coding(gff, [("chr", 120)], strand) == WRAPPED_CDS
        gff = tmp_path / f"noregion{strand}.gff3"
        gff.write_text("##gff-version 3\n" + "".join(lines))
        assert _wrapped_coding(gff, [("chr", 120)], strand) == WRAPPED_CDS  # The reference's length
        # Neither (the annotated sequence is taken for the reference by its length... which is unknown): the parts
        # are read in coordinate order, as documented
        assert _wrapped_coding(gff, [("other", 120)], strand) == WRAPPED_CDS[30:] + WRAPPED_CDS[:30]
    # Three parts: the groups either side of the origin are split at the widest gap
    gff = tmp_path / "three.gff3"
    gff.write_text("##gff-version 3\n##sequence-region chr 1 120\n"
                   "chr\t.\tCDS\t1\t30\t.\t+\t0\tID=c1\nchr\t.\tCDS\t80\t84\t.\t+\t0\tID=c1\n"
                   "chr\t.\tCDS\t91\t120\t.\t+\t0\tID=c1\n")
    ann = load_annotation(gff, [("chr", 120)])
    assert ann.sequences["chr"].genes[0].cds[0].parts == [(80, 84), (91, 120), (1, 30)]


def test_gff3_single_line_feature_across_the_origin(tmp_path):
    # Bakta writes a feature across the origin of a circular contig as one line whose end is beyond the length
    # (end = its end + the length): split into its two parts, in the order of translation
    for strand in "+-":
        gff = tmp_path / f"bakta{strand}.gff3"
        gff.write_text("##gff-version 3\n##sequence-region chr 1 120\n"
                       "chr\tBakta\tregion\t1\t120\t.\t+\t.\tID=chr;Name=chr;Is_circular=true\n"
                       f"chr\tBakta\tCDS\t91\t150\t.\t{strand}\t0\tID=c1;gene=wrap\n")
        assert _wrapped_coding(gff, [("chr", 120)], strand) == WRAPPED_CDS
    gff = tmp_path / "linear.gff3"  # Not flagged circular: beyond the end of the sequence, as before
    gff.write_text("##gff-version 3\n##sequence-region chr 1 120\nchr\tx\tCDS\t91\t150\t.\t+\t0\tID=c1;gene=wrap\n")
    ann = load_annotation(gff, [("chr", 120)])
    assert not ann.sequences and any("beyond the end" in w for w in ann.warnings)


# ---------------------------------------------------------------------------------------------------------------
# Gene labels: a gene without a symbol is shown as "locus_tag (product)"
# ---------------------------------------------------------------------------------------------------------------

def test_gene_without_symbol_is_labelled_by_locus_tag_and_product():
    rec = read_genbank_text(_gb("""
    gene            1..30
                    /locus_tag="LK299_pgp087"
    CDS             1..30
                    /locus_tag="LK299_pgp087"
                    /product="maturase K"
    gene            complement(41..60)
                    /gene="trnX"
                    /locus_tag="SIM_t001"
    tRNA            complement(41..60)
                    /gene="trnX"
                    /product="tRNA-Xaa"
    gene            70..84
                    /locus_tag="SIM_p010"
    CDS             70..84
                    /locus_tag="SIM_p010"
                    /product="ribosomal protein S12"
    CDS             90..110
                    /product="hypothetical protein"
    gene            111..120
                    /locus_tag="SIM_p011"
    CDS             111..120
                    /gene="ycf1"
                    /locus_tag="SIM_p011"
                    /product="Ycf1"
    """))
    genes = {g.name: g for g in genes_from_genbank(rec.features)}
    unnamed = genes["LK299_pgp087"]  # The identifier is still the locus tag (nothing keyed on it changes)
    assert unnamed.symbol == "" and unnamed.product == "maturase K"
    assert unnamed.label == "LK299_pgp087 (maturase K)" and unnamed.map_label == "maturase K"  # 10 characters
    named = genes["trnX"]
    assert named.symbol == "trnX" and named.label == "trnX" and named.map_label == "trnX"
    long_product = genes["SIM_p010"]
    assert long_product.label == "SIM_p010 (ribosomal protein S12)" and long_product.map_label == "SIM_p010"
    only_product = genes["hypothetical protein"]  # Named after its product: not repeated
    assert only_product.label == "hypothetical protein" and only_product.map_label == "hypothetical protein"
    cds_symbol = genes["SIM_p011"]  # The gene feature has no /gene but its CDS has: the symbol is shown
    assert cds_symbol.symbol == "ycf1" and cds_symbol.label == "SIM_p011 (ycf1)"
    assert cds_symbol.map_label == "SIM_p011"
    # The intergenic context names the neighbours by their labels
    from bacon.annotation import SequenceAnnotation, annotate_snp
    seq_ann = SequenceAnnotation("ref", 120, list(genes.values()), [], True)
    assert annotate_snp(seq_ann, 35, "A", ["C"], SEQ).context == \
        "intergenic between LK299_pgp087 (maturase K) and trnX"


def test_gff3_symbols_are_gene_attributes_or_names_that_are_not_identifiers(tmp_path):
    gff = tmp_path / "symbols.gff3"
    gff.write_text(
        "##gff-version 3\n##sequence-region chr 1 1000\n"
        # NCBI, no symbol: Name is the locus tag; the CDS's Name is its protein accession
        "chr\t.\tgene\t1\t30\t.\t+\t.\tID=gene-LK299_pgp087;Name=LK299_pgp087;locus_tag=LK299_pgp087\n"
        "chr\t.\tCDS\t1\t30\t.\t+\t0\tID=cds-YP_1;Parent=gene-LK299_pgp087;Name=YP_1;locus_tag=LK299_pgp087;"
        "product=maturase K\n"
        # NCBI, a symbol
        "chr\t.\tgene\t41\t60\t.\t+\t.\tID=gene-SIM_p010;Name=matK;gene=matK;locus_tag=SIM_p010\n"
        "chr\t.\tCDS\t41\t60\t.\t+\t0\tID=cds-YP_2;Parent=gene-SIM_p010;gene=matK;product=maturase K\n"
        # A Name that is a symbol (no gene=, no locus tag), as Ensembl writes one
        "chr\t.\tgene\t101\t130\t.\t+\t.\tID=gene:AT1G01010;Name=NAC001;biotype=protein_coding\n"
        "chr\t.\tCDS\t101\t130\t.\t+\t0\tID=CDS:AT1G01010.1;Parent=gene:AT1G01010;product=NAC domain protein\n"
        # A Name that is the ID without its type prefix: a symbol unless the gene has a locus tag or the Name looks
        # like one (PREFIX_number); a Name that is the product: not a symbol
        "chr\t.\tgene\t201\t230\t.\t+\t.\tID=gene-G3;Name=G3\n"
        "chr\t.\tCDS\t201\t230\t.\t+\t0\tID=cds-G3;Parent=gene-G3;product=hypothetical protein\n"
        "chr\t.\tgene\t231\t260\t.\t+\t.\tID=gene-matK;Name=matK\n"
        "chr\t.\tCDS\t231\t260\t.\t+\t0\tID=cds-matK;Parent=gene-matK;product=maturase K\n"
        "chr\t.\tgene\t261\t290\t.\t+\t.\tID=gene-SIM_p020;Name=SIM_p020\n"
        "chr\t.\tCDS\t261\t290\t.\t+\t0\tID=cds-SIM_p020;Parent=gene-SIM_p020;product=hypothetical protein\n"
        "chr\t.\tgene\t291\t299\t.\t+\t.\tID=gene-G5;Name=G5;locus_tag=L_G5\n"
        "chr\t.\tCDS\t291\t299\t.\t+\t0\tID=cds-G5;Parent=gene-G5;product=hypothetical protein\n"
        "chr\t.\tgene\t301\t330\t.\t+\t.\tID=g4;Name=hypothetical protein;product=hypothetical protein\n"
        # A CDS without a gene feature: its Name is not a symbol either
        "chr\t.\tCDS\t401\t430\t.\t+\t0\tID=cds-YP_5;Name=YP_5;locus_tag=L5;product=photosystem I protein\n")
    genes = {g.name: g for g in genes_from_gff3(read_gff3(gff)[0])}
    assert genes["LK299_pgp087"].symbol == "" and genes["LK299_pgp087"].label == "LK299_pgp087 (maturase K)"
    assert genes["matK"].symbol == "matK" and genes["matK"].label == "matK"
    assert genes["NAC001"].symbol == "NAC001" and genes["NAC001"].label == "NAC001"
    assert genes["G3"].symbol == "G3" and genes["G3"].label == "G3"  # No locus tag, not shaped like one
    assert genes["matK"].symbol == "matK" and genes["matK"].label == "matK"
    assert genes["SIM_p020"].symbol == "" and genes["SIM_p020"].label == "SIM_p020 (hypothetical protein)"
    assert genes["L_G5"].symbol == "" and genes["L_G5"].label == "L_G5 (hypothetical protein)"  # Has a locus tag
    assert genes["hypothetical protein"].symbol == "" and genes["hypothetical protein"].label == "hypothetical protein"
    assert genes["L5"].symbol == "" and genes["L5"].label == "L5 (photosystem I protein)"


# ---------------------------------------------------------------------------------------------------------------
# The inverted repeat detected in the sequence
# ---------------------------------------------------------------------------------------------------------------

def make_plastome(seed: int = 1, lsc: int = 9000, ir: int = 6000, ssc: int = 4000, start: int = 0,
                  mismatches: int = 0, deletion: int = 0, duplication: bool = False, spacing: int = 7,
                  edit=None) -> str:
    """A circular plastome-like random sequence LSC + IRb + SSC + IRa (IRa the reverse complement of IRb), whose
    junction bases do not match across (so the repeat ends exactly where it is planted), rotated to start at
    `start` (0-based). `mismatches` bases of IRa are changed (from its base 100, every `spacing` bases), and
    `deletion` bases removed from its middle; with `duplication`, the repeat contains a 40 bp stretch duplicated
    50 bp further (in both copies); `edit`, a function, rewrites IRa last."""
    import random

    from bacon.annotation import reverse_complement
    rng = random.Random(seed)
    parts = ["".join(rng.choices("ACGT", k=n)) for n in (lsc, ir, ssc)]
    lsc_seq, irb, ssc_seq = parts
    lsc_seq = "A" + lsc_seq[1:-1] + "A"  # Across the origin: A pairs with T, not A
    ssc_seq = "A" + ssc_seq[1:-1] + "A"
    if duplication:
        irb = irb[:1050] + irb[1000:1040] + irb[1090:]
    ira = list(reverse_complement(irb))
    for p in range(100, 100 + spacing * mismatches, spacing):
        ira[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[ira[p]]
    if deletion:
        del ira[ir // 2:ir // 2 + deletion]
    ira = "".join(ira)
    seq = lsc_seq + irb + ssc_seq + (edit(ira) if edit else ira)
    return seq[start:] + seq[:start]


def test_detects_the_inverted_repeat_and_derives_the_regions():
    from bacon.annotation import detect_inverted_repeat, regions_from_repeat
    seq = make_plastome()
    repeat = detect_inverted_repeat(seq)
    assert repeat is not None
    assert repeat.first == (9001, 15000) and repeat.second == (19001, 25000)
    assert repeat.lengths == (6000, 6000) and repeat.differences == 0 and repeat.identity == 1.0
    assert repeat.text() == "two copies of 6,000 bp, 100% identical"
    assert [(r.name, r.start, r.end, r.length) for r in regions_from_repeat(repeat, len(seq))] == [
        ("LSC", 1, 9000, 9000), ("IRb", 9001, 15000, 6000), ("SSC", 15001, 19000, 4000), ("IRa", 19001, 25000, 6000)]
    assert detect_inverted_repeat(seq.lower()) == repeat  # Case does not matter


def test_detection_tolerates_mismatches_and_a_small_deletion():
    from bacon.annotation import detect_inverted_repeat
    repeat = detect_inverted_repeat(make_plastome(mismatches=5))
    assert repeat is not None and repeat.first == (9001, 15000) and repeat.second == (19001, 25000)
    assert repeat.differences == 5 and 0.999 < repeat.identity < 1
    assert repeat.text() == "two copies of 6,000 bp, 99.92% identical"
    repeat = detect_inverted_repeat(make_plastome(deletion=3, mismatches=1))
    assert repeat is not None and repeat.first == (9001, 15000) and repeat.second == (19001, 24997)
    assert repeat.lengths == (6000, 5997) and repeat.differences == 2  # The deletion is one difference
    assert repeat.text() == "copies of 6,000 and 5,997 bp, 99.97% identical"
    # Too many differences along the whole repeat: not the near-identical repeat of a plastome
    assert detect_inverted_repeat(make_plastome(mismatches=70, spacing=80)) is None  # 70 / 6000 > 1%
    # The same 70 mismatches packed into 490 bp near one end: a diverged flank, dropped (the chain is trimmed to
    # its part scoring most); the 5.4 kb left are a repeat (the Arabidopsis mitochondrion NC_037304.1 has a
    # 6,590 bp repeat, 100% identical, in diverged flanks)
    repeat = detect_inverted_repeat(make_plastome(mismatches=70))
    assert repeat is not None and repeat.first == (9001, 14416) and repeat.second == (19585, 25000)
    assert repeat.lengths == (5416, 5416) and repeat.differences == 0
    # A short duplication inside the repeat seeds a parallel antidiagonal (a k-mer matching 50 bp further in the
    # other copy); it is not read as two 50 bp indels (the rice NC_001320.1 repeat has one)
    repeat = detect_inverted_repeat(make_plastome(duplication=True))
    assert repeat is not None and repeat.lengths == (6000, 6000) and repeat.differences == 0


def test_detection_across_the_origin_and_regions_naming():
    from bacon.annotation import detect_inverted_repeat, regions_from_repeat
    n = 25000
    for start in (22000, 9000 + 5999, 12000):  # IRa across the origin; IRb's last base first; starting in IRb
        seq = make_plastome(start=start)

        def at(pos: int, start: int = start) -> int:
            return (pos - 1 - start) % n + 1

        repeat = detect_inverted_repeat(seq)
        assert repeat is not None, start
        copies = sorted((repeat.first, repeat.second))
        assert copies == sorted(((at(9001), at(15000)), (at(19001), at(25000)))), start
        assert repeat.lengths == (6000, 6000) and repeat.differences == 0
        regions = {r.name: (r.start, r.end, r.length) for r in regions_from_repeat(repeat, n)}
        assert regions == {"LSC": (at(1), at(9000), 9000), "IRb": (at(9001), at(15000), 6000),
                           "SSC": (at(15001), at(19000), 4000), "IRa": (at(19001), at(25000), 6000)}, start
    # A deletion and a mismatch in the part of IRa across the origin: the extension through the origin stops at
    # the deletion, and the search again from between the copies gets the whole repeat
    repeat = detect_inverted_repeat(make_plastome(start=22000, deletion=3, mismatches=1))
    assert repeat is not None and repeat.lengths == (6000, 5997) and repeat.differences == 2
    assert (repeat.first, repeat.second) == ((11998, 17997), (21998, 2997))  # IRb, then IRa across the origin
    # Starting in the SSC: the gap between the copies is the LSC, so the first copy is IRa
    seq = make_plastome(start=17000)
    names = [r.name for r in regions_from_repeat(detect_inverted_repeat(seq), n)]
    assert names == ["IRa", "LSC", "IRb", "SSC"]


def test_no_repeat_detected_without_one():
    import random
    import time

    from bacon.annotation import MAX_DETECTION_LENGTH, detect_inverted_repeat
    rng = random.Random(7)
    assert detect_inverted_repeat("".join(rng.choices("ACGT", k=25000))) is None  # No repeat
    assert detect_inverted_repeat(make_plastome(ir=3000)) is None  # Below the minimum (the bundled example's size)
    assert detect_inverted_repeat(make_plastome(ir=3000, lsc=2000, ssc=1000)) is None  # Too short to search
    assert detect_inverted_repeat("ACGT" * 10000) is None  # Low complexity: its own reverse complement everywhere
    half = "".join(rng.choices("ACGT", k=6000))
    from bacon.annotation import find_regions, regions_from_repeat, reverse_complement
    palindrome = "".join(rng.choices("ACGT", k=8000)) + half + reverse_complement(half) \
        + "".join(rng.choices("ACGT", k=8000))
    # The two halves of a palindrome are two copies that abut (the Toxoplasma apicoplast NC_001799.1 has its
    # copies abutting across the origin): a repeat, but no regions, and a band without regions saying why
    repeat = detect_inverted_repeat(palindrome)
    assert repeat is not None and (repeat.first, repeat.second) == ((8001, 14000), (14001, 20000))
    assert regions_from_repeat(repeat, len(palindrome)) == []
    band = find_regions([], len(palindrome), palindrome)
    assert band.regions == [] and band.source == "none" and band.repeat == repeat
    assert band.note == "inverted repeat found but not a plastome layout (the copies abut)"
    assert band.record() == {"source": "none", "regions": [], "note": band.note,
                             "repeat": {"copies": [[8001, 14000], [14001, 20000]], "lengths": [6000, 6000],
                                        "differences": 0, "identity": 1.0}}
    assert detect_inverted_repeat("A" * (MAX_DETECTION_LENGTH + 1)) is None  # Not searched
    bacterial = "".join(rng.choices("ACGT", k=1_000_000))
    started = time.perf_counter()
    assert detect_inverted_repeat(bacterial) is None
    assert time.perf_counter() - started < 5  # About 0.4 s; a plastome takes a tenth of a second
    # N bases are not matched: a repeat of N is not a repeat
    assert detect_inverted_repeat("".join(rng.choices("ACGT", k=8000)) + "N" * 6000
                                  + "".join(rng.choices("ACGT", k=4000)) + "N" * 6000) is None


def test_find_regions_prefers_the_annotation_unless_it_contradicts_the_sequence():
    from bacon.annotation import find_regions
    seq = make_plastome()
    n = len(seq)
    annotated = _features("""
    repeat_region   9001..15000
                    /rpt_type=inverted
    repeat_region   19001..25000
                    /rpt_type=inverted
    """)
    warnings: list[str] = []
    band = find_regions(annotated, n, seq, "ref", warnings)
    assert band.source == "annotation" and band.repeat is None and not warnings
    assert band.text() == "the annotated inverted repeats"
    # Annotated a little differently (within 20%): still the annotation's coordinates
    shifted = _features("""
    repeat_region   9201..15000
                    /note="IRb"
    repeat_region   19001..24800
                    /note="IRa"
    """)
    band = find_regions(shifted, n, seq, "ref", warnings)
    assert band.source == "annotation" and [r.start for r in band.regions] == [9201, 15001, 19001, 24801]
    # Other features named as the repeats (NC_001879.2 annotates its LSC as "inverted repeat B"): the sequence wins
    wrong = _features("""
    repeat_region   1..9000
                    /note="inverted repeat B"
    repeat_region   15001..19000
                    /note="inverted repeat A"
    """)
    band = find_regions(wrong, n, seq, "ref", warnings)
    assert band.source == "sequence" and [r.start for r in band.regions] == [1, 9001, 15001, 19001]
    assert warnings == ["ref: the annotated inverted repeats (1–9,000, 15,001–19,000) are not the inverted repeat "
                        "found in the sequence (9,001–15,000 and 19,001–25,000); the regions follow the sequence"]
    # No annotated repeats: the sequence; no sequence: the annotation only; neither: nothing
    assert find_regions([], n, seq).source == "sequence"
    assert find_regions(wrong, n, None).source == "annotation"
    assert find_regions([], n, None) is None and find_regions([], 25000, "ACGT" * 6250) is None
    record = find_regions([], n, seq).record()
    assert record["source"] == "sequence" and record["regions"][1] == {"name": "IRb", "start": 9001, "end": 15000,
                                                                       "length": 6000}
    assert record["repeat"] == {"copies": [[9001, 15000], [19001, 25000]], "lengths": [6000, 6000],
                                "differences": 0, "identity": 1.0}


def test_load_annotation_detects_the_regions_when_given_the_sequences(tmp_path):
    seq = make_plastome()
    gb = tmp_path / "plastid.gb"
    gb.write_text(_gb("""
    gene            101..400
                    /gene="psbA"
    CDS             101..400
                    /gene="psbA"
    """, seq=seq, name="plastid", version="plastid.1"))
    without = load_annotation(gb, [("plastid.1", len(seq))])
    assert without.sequences["plastid.1"].regions == [] and not without.has_regions
    with_seq = load_annotation(gb, [("plastid.1", len(seq))], seqs={"plastid.1": seq})
    band = with_seq.sequences["plastid.1"].band
    assert with_seq.has_regions and band.source == "sequence" and band.repeat.lengths == (6000, 6000)
    assert [r.name for r in with_seq.sequences["plastid.1"].regions] == ["LSC", "IRb", "SSC", "IRa"]
    assert with_seq.sequences["plastid.1"].region_at(12000) == "IRb"


# ---------------------------------------------------------------------------------------------------------------
# Detection: the differences between the copies, indels and N runs, flanks, abutting copies, the rotated search
# ---------------------------------------------------------------------------------------------------------------

def _random(n: int, seed: int) -> str:
    import random
    return "".join(random.Random(seed).choices("ACGT", k=n))


def _edit_distance(a: str, b: str) -> int:
    """Plain edit distance (reference for the banded one)."""
    previous = list(range(len(b) + 1))
    for x, ca in enumerate(a, 1):
        current = [x]
        for y, cb in enumerate(b, 1):
            current.append(min(previous[y] + 1, current[y - 1] + 1, previous[y - 1] + (ca != cb)))
        previous = current
    return previous[-1]


def test_banded_edit_distance_matches_the_plain_one_within_the_band():
    import random

    from bacon.annotation import _banded_edit_distance
    rng = random.Random(3)
    for _ in range(60):
        a = "".join(rng.choices("ACGT", k=rng.randint(1, 60)))
        b = list(a)
        for _ in range(rng.randint(0, 6)):  # A few substitutions and indels (within a band of 8)
            pos = rng.randrange(len(b) + 1)
            op = rng.choice("sid")
            if op == "s" and pos < len(b):
                b[pos] = rng.choice("ACGT")
            elif op == "i":
                b.insert(pos, rng.choice("ACGT"))
            elif op == "d" and pos < len(b) and len(b) > 1:
                del b[pos]
        b = "".join(b)
        short, long_ = (a, b) if len(a) <= len(b) else (b, a)
        assert _banded_edit_distance(short, long_, 8) == _edit_distance(a, b), (a, b)
    assert _banded_edit_distance("", "ACGT", 2) == 4 and _banded_edit_distance("ACGT", "ACGT", 0) == 0
    assert _banded_edit_distance("ACNT", "ACGTT", 2) == 1  # N matches anything


def test_stretch_differences_count_events_not_bases():
    import bacon.annotation as module
    from bacon.annotation import _stretch_differences
    x = _random(300, 5)
    assert _stretch_differences(x, x) == 0 and _stretch_differences("", "") == 0
    assert _stretch_differences(x, x[1:] + "G") == 2  # A deletion and an insertion 300 bp apart: not 225 mismatches
    assert _stretch_differences(x, x[:100] + "GGGGG" + x[100:]) == 1  # An insertion of 5 bases: one difference
    assert _stretch_differences(x, "") == 1 and _stretch_differences("", "ACGT") == 1
    two = x[:50] + ("A" if x[50] != "A" else "C") + x[51:200] + ("A" if x[200] != "A" else "C") + x[201:]
    assert _stretch_differences(x, two) == 2  # Mismatches, compared base by base
    ten = list(x)
    for p in range(10, 300, 29):
        ten[p] = "A" if ten[p] != "A" else "C"
    assert _stretch_differences(x, "".join(ten)) == 10  # More than IR_FEW: aligned, still 10
    assert _stretch_differences(x, x[:100] + "N" * 50 + x[150:]) == 1  # A run of N: one difference
    assert _stretch_differences(x, x[:100] + "N" * 50 + x[100:]) == 2  # Inserted: an indel and a run of N
    r = x[:10] + "R" + x[11:]
    assert _stretch_differences(r, r) == 0 and _stretch_differences(x, r) == 1  # The same code in both: nothing
    # A long insertion between two mismatches: aligned within a band of its length (plus the margin); an
    # alignment too large to do counts the stretch as all different (one more than its bases)
    big = x[:150] + "T" * 1500 + x[150:]
    assert _stretch_differences(x, big) == 1  # Only the insertion is left once the common ends are removed
    assert _stretch_differences(two, big) == 3
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(module, "IR_MAX_CELLS", 10)
        assert _stretch_differences(two, big) == 152  # The 151 bases from one mismatch to the other, plus one


def test_compensating_indels_do_not_reject_the_repeat():
    # A base deleted and one inserted 300 bp further in IRa (NC_007144.1, cucumber, has such a pair): the stretch
    # between, of equal length in both copies, differs at three bases in four read base by base
    from bacon.annotation import detect_inverted_repeat
    repeat = detect_inverted_repeat(make_plastome(edit=lambda ira: ira[:3000] + ira[3001:3300] + "G" + ira[3300:]))
    assert repeat is not None and repeat.first == (9001, 15000) and repeat.second == (19001, 25000)
    assert repeat.lengths == (6000, 6000) and repeat.differences == 2


def test_insertions_and_n_runs_in_one_copy():
    from bacon.annotation import IR_MAX_INDEL, detect_inverted_repeat
    for size in (100, 600, 1500, 2500):  # Larger than the band of 500 antidiagonals: the chain goes on past it
        repeat = detect_inverted_repeat(make_plastome(edit=lambda ira, s=size: ira[:3000] + _random(s, 9) + ira[3000:]))
        assert repeat is not None and repeat.first == (9001, 15000) and repeat.second == (19001, 25000 + size)
        assert repeat.lengths == (6000, 6000 + size) and repeat.differences == 1, size  # One difference
    assert IR_MAX_INDEL == 3000
    # An insertion larger than IR_MAX_INDEL cuts the repeat: the longer piece is the repeat
    seq = make_plastome(lsc=12000, ir=12000, edit=lambda ira: ira[:5000] + _random(3500, 9) + ira[5000:])
    repeat = detect_inverted_repeat(seq)
    assert repeat is not None and repeat.lengths == (7000, 7000) and repeat.differences == 0
    assert repeat.first == (12001, 19000) and repeat.second == (36501, 43500)
    # A run of 300 N in one copy (a scaffold gap) is one difference; a run of 150 at the very start of IRa too
    repeat = detect_inverted_repeat(make_plastome(edit=lambda ira: ira[:3000] + "N" * 300 + ira[3300:]))
    assert repeat is not None and repeat.lengths == (6000, 6000) and repeat.differences == 1
    repeat = detect_inverted_repeat(make_plastome(edit=lambda ira: ira[:3000] + "N" * 300 + ira[3000:]))
    assert repeat is not None and repeat.lengths == (6000, 6300) and repeat.differences == 1  # Inserted: an indel
    # A run of N longer than IR_MAX_GAP leaves no seed on either side of it: the repeat is cut there
    seq = make_plastome(lsc=12000, ir=12000, edit=lambda ira: ira[:4000] + "N" * 2500 + ira[6500:])
    repeat = detect_inverted_repeat(seq)
    assert repeat is not None and repeat.lengths == (5500, 5500) and repeat.first == (12001, 17500)


def test_a_truncated_detection_confirms_a_correct_annotation():
    # The annotation gives the whole copies; the detection stops at an insertion larger than IR_MAX_INDEL (or
    # at a run of N longer than IR_MAX_GAP): each detected copy lies inside an annotated one, so the annotation
    # is kept. Only copies found mostly outside the annotated ones override them.
    from bacon.annotation import _agree, find_regions
    seq = make_plastome(lsc=12000, ir=12000, edit=lambda ira: ira[:5000] + _random(3500, 9) + ira[5000:])
    annotated = _features("""
    repeat_region   12001..24000
                    /rpt_type=inverted
    repeat_region   28001..43500
                    /rpt_type=inverted
    """)
    warnings: list[str] = []
    band = find_regions(annotated, len(seq), seq, "ref", warnings)
    assert band.source == "annotation" and not warnings
    assert [(r.name, r.start, r.end) for r in band.regions] == [("LSC", 1, 12000), ("IRb", 12001, 24000),
                                                                ("SSC", 24001, 28000), ("IRa", 28001, 43500)]
    seq = make_plastome(lsc=12000, ir=12000, edit=lambda ira: ira[:4000] + "N" * 2500 + ira[6500:])
    annotated = _features("""
    repeat_region   12001..24000
                    /rpt_type=inverted
    repeat_region   28001..40000
                    /rpt_type=inverted
    """)
    band = find_regions(annotated, len(seq), seq, "ref", warnings)
    assert band.source == "annotation" and not warnings
    # The threshold: a detected copy lying 80% inside an annotated copy confirms it, 78% does not
    from bacon.annotation import InvertedRepeat, Region
    regions = [Region("IRb", 9001, 15000, 6000), Region("IRa", 19001, 25000, 6000)]
    assert _agree(regions, InvertedRepeat((10201, 16200), (17801, 23800), (6000, 6000), 0), 25000)  # 4,800 inside
    assert not _agree(regions, InvertedRepeat((10301, 16300), (17801, 23800), (6000, 6000), 0), 25000)  # 4,700
    assert not _agree(regions, InvertedRepeat((10201, 16200), (13801, 19800), (6000, 6000), 0), 25000)  # One copy


def test_abutting_copies_across_the_origin_and_the_rotated_search():
    # The Toxoplasma apicoplast NC_001799.1 has its copies abutting across the origin: a repeat with no region
    # between the copies on one side, so no band; recorded with a note
    import bacon.annotation as module
    from bacon.annotation import _same_pair, detect_inverted_repeat, find_regions, reverse_complement
    x = _random(5500, 11)
    seq = reverse_complement(x) + "A" + _random(19998, 12) + "A" + x  # A pairs with T: the copies end there
    repeat = detect_inverted_repeat(seq)
    assert repeat is not None and (repeat.first, repeat.second) == ((1, 5500), (25501, 31000))
    band = find_regions([], len(seq), seq)
    assert band.regions == [] and band.source == "none" and band.repeat == repeat
    assert band.note == "inverted repeat found but not a plastome layout (the copies abut)"
    # The search in the rotated sequence (a second, equal search) only happens when a copy crosses the origin, or
    # ends at it with the repeat going on beyond: not for the common layout where IRa ends at the last base
    calls: list[int] = []
    seeds = module._repeat_seeds
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(module, "_repeat_seeds", lambda s, k, step: calls.append(len(s)) or seeds(s, k, step))
        assert detect_inverted_repeat(make_plastome()).second == (19001, 25000) and len(calls) == 1
        calls.clear()
        assert detect_inverted_repeat(make_plastome(start=22000)) is not None and len(calls) == 2  # Across
        calls.clear()
        assert detect_inverted_repeat(make_plastome(start=9000 + 5999)).lengths == (6000, 6000) and len(calls) == 2
        calls.clear()
        assert detect_inverted_repeat(seq) == repeat and len(calls) == 1  # Abutting: no room to go on
    # The rotated search must find the same pair of copies (compared in either pairing: a copy across the origin
    # starts late in the sequence but covers its first base)
    from bacon.annotation import InvertedRepeat
    a = InvertedRepeat((1, 2997), (12001, 14997), (2997, 2997), 0)
    b = InvertedRepeat((12001, 18000), (22001, 2997), (6000, 5997), 2)
    assert _same_pair(b, a, 24997) and _same_pair(a, b, 24997)
    assert not _same_pair(InvertedRepeat((100, 6000), (30000, 36000), (5901, 6001), 0), b, 40000)


def test_identity_is_recorded_rounded():
    from bacon.annotation import find_regions
    seq = make_plastome(mismatches=5)
    record = find_regions([], len(seq), seq).record()
    assert record["repeat"]["differences"] == 5 and record["repeat"]["identity"] == 0.999167  # Not 0.99916666...


# ---------------------------------------------------------------------------------------------------------------
# A plastome layout, or not
# ---------------------------------------------------------------------------------------------------------------

def test_a_repeat_that_is_not_a_plastome_layout_gives_no_band():
    from bacon.annotation import PLASTOME_MAX_SINGLE_COPY, PLASTOME_MIN_IR_FRACTION, find_regions
    assert PLASTOME_MIN_IR_FRACTION == 0.05 and PLASTOME_MAX_SINGLE_COPY == 200_000
    # Two copies of 5 kb in 205 kb (4.9%): the inverted rRNA operons of a bacterium (H. pylori NC_000915.1 has
    # 10.5 kb copies in 1.67 Mb) rather than the repeats of a plastome
    seq = make_plastome(lsc=190000, ir=5000, ssc=5000)
    band = find_regions([], len(seq), seq)
    assert band.source == "none" and band.regions == [] and band.repeat.lengths == (5000, 5000)
    assert band.note == "inverted repeat found but not a plastome layout (the repeats are 4.9% of the sequence, " \
                        "less than 5%)"
    # A single-copy region of 205 kb (the Zea mays mitochondrion NC_007982.1 has 481 kb between 16.9 kb copies)
    seq = make_plastome(lsc=205000, ir=12000, ssc=5000)
    band = find_regions([], len(seq), seq)
    assert band.source == "none" and band.note == "inverted repeat found but not a plastome layout (the larger " \
                                                  "single-copy region is 205 kb, more than 200 kb)"
    # Just plausible: 5% and 200 kb
    seq = make_plastome(lsc=200000, ir=6000, ssc=28000)
    assert [r.name for r in find_regions([], len(seq), seq).regions] == ["LSC", "IRb", "SSC", "IRa"]


def test_a_sequence_said_not_to_be_a_plastid_gets_no_plastome_regions():
    from bacon.annotation import Location, RawFeature, find_regions
    seq = make_plastome()
    n = len(seq)
    mito = _features(f"""
    source          1..{n}
                    /organelle="mitochondrion"
    repeat_region   9001..15000
                    /rpt_type=inverted
    repeat_region   19001..25000
                    /rpt_type=inverted
    """)
    warnings: list[str] = []
    band = find_regions(mito, n, seq, "mt", warnings)  # The Arabidopsis mitochondrion NC_037304.1 annotates a pair
    assert band.source == "none" and band.regions == [] and band.repeat.lengths == (6000, 6000)
    assert band.note == "inverted repeat found but not a plastome layout (the source feature says " \
                        "organelle=mitochondrion)"
    assert warnings == ["mt: the annotated inverted repeats (9,001–15,000, 19,001–25,000) are not a plastome layout "
                        "(the source feature says organelle=mitochondrion); no regions from them"]
    band = find_regions(mito, n, None, "mt", warnings)  # Without the sequence: the annotation alone, rejected
    assert band.source == "none" and band.repeat is None and band.regions == []
    assert band.note == "annotated inverted repeats but not a plastome layout (the source feature says " \
                        "organelle=mitochondrion)"
    assert band.record() == {"source": "none", "regions": [], "note": band.note}
    plastid = _features(f"""
    source          1..{n}
                    /organelle="plastid:chloroplast"
    """)
    assert find_regions(plastid, n, seq).source == "sequence"
    # GFF3: NCBI's region feature says genome=chromosome (a bacterium), chloroplast, mitochondrion
    for genome, expected in (("chromosome", "none"), ("chloroplast", "sequence"), ("mitochondrion", "none")):
        region = RawFeature("region", "ref", Location(1, [(1, n)]), {"genome": genome, "mol_type": "genomic DNA"})
        band = find_regions([region], n, seq)
        assert band.source == expected, genome
        if expected == "none":
            assert band.note.endswith(f"(the region feature says genome={genome})")


def test_annotated_repeats_that_give_no_regions_are_said_so():
    from bacon.annotation import find_regions
    seq = make_plastome()
    n = len(seq)
    # NC_001879.2 (tobacco) annotates its LSC as `inverted repeat B`, its SSC as `SSC; inverted repeat A` and its
    # IRa as `IRA`: the two IRa features merge into one touching the IRb feature at the origin, so no regions
    tobacco = _features("""
    repeat_region   1..9000
                    /note="inverted repeat B"
    repeat_region   15001..19000
                    /note="SSC; inverted repeat A"
    misc_feature    19001..25000
                    /note="IRA"
    """)
    assert derive_regions(tobacco, n) == []
    warnings: list[str] = []
    band = find_regions(tobacco, n, seq, "NC_001879.2", warnings)
    assert band.source == "sequence" and [r.start for r in band.regions] == [1, 9001, 15001, 19001]
    assert warnings == ["NC_001879.2: the annotated inverted repeats do not give a plastome layout; the regions "
                        "follow the inverted repeat found in the sequence (9,001–15,000 and 19,001–25,000)"]
    assert find_regions(tobacco, n, None, "x", warnings) is None and len(warnings) == 1  # Nothing to say
    # Annotated repeats too small for a plastome layout (two of 600 bp: 4.8% of the sequence), the sequence has one
    tiny = _features("""
    repeat_region   9001..9600
                    /rpt_type=inverted
    repeat_region   24401..25000
                    /rpt_type=inverted
    """)
    warnings.clear()
    band = find_regions(tiny, n, seq, "ref", warnings)
    assert band.source == "sequence" and band.repeat.lengths == (6000, 6000)
    assert warnings == ["ref: the annotated inverted repeats (9,001–9,600, 24,401–25,000) are not a plastome layout "
                        "(the repeats are 4.8% of the sequence, less than 5%); the regions follow the inverted "
                        "repeat found in the sequence (9,001–15,000 and 19,001–25,000)"]


# ---------------------------------------------------------------------------------------------------------------
# Reading: wrapped GenBank values, Ensembl and NCBI GFF3 oddities
# ---------------------------------------------------------------------------------------------------------------

def test_genbank_wrapped_values_join_without_a_space_only_inside_a_broken_word():
    rec = read_genbank_text(_gb("""
    CDS             1..30
                    /product="2,
                    6-diaminopimelate ATP-
                    dependent protein A,
                    chloroplastic"
                    /note="some words
                    more words, 3 of them"
    """))
    q = rec.features[0].qualifiers
    assert q["product"] == "2,6-diaminopimelate ATP-dependent protein A, chloroplastic"
    assert q["note"] == "some words more words, 3 of them"


def test_gff3_ensembl_ncrna_genes_and_ncbi_parentless_rnas(tmp_path):
    gff = tmp_path / "odd.gff3"
    gff.write_text(
        "##gff-version 3\n##sequence-region Pt 1 1000\n"
        # Ensembl: a tRNA gene is an ncRNA_gene, without a Name when it has no symbol; its gene_id names it
        "Pt\tara\tncRNA_gene\t4\t76\t.\t-\t.\tID=gene:ATCG00010;biotype=tRNA;gene_id=ATCG00010\n"
        "Pt\tara\ttRNA\t4\t76\t.\t-\t.\tID=transcript:ATCG00010.1;Parent=gene:ATCG00010;biotype=tRNA\n"
        "Pt\tara\texon\t4\t76\t.\t-\t.\tParent=transcript:ATCG00010.1;Name=ATCG00010.1.exon1\n"
        "Pt\tara\tncRNA_gene\t101\t200\t.\t+\t.\tID=gene:ATCG00920;Name=RRN16S;biotype=rRNA;gene_id=ATCG00920\n"
        "Pt\tara\trRNA\t101\t200\t.\t+\t.\tID=transcript:ATCG00920.1;Parent=gene:ATCG00920;biotype=rRNA\n"
        # NCBI: a gene with gene_biotype=other whose rRNA has no Parent (tomato NC_007898.3), and a parentless
        # rRNA with another locus tag, and one on the other strand: genes of their own
        "Pt\tRefSeq\tgene\t301\t420\t.\t-\t.\tID=gene-LyesC2r005;Name=LyesC2r005;gene_biotype=other;"
        "locus_tag=LyesC2r005\n"
        "Pt\tRefSeq\trRNA\t301\t420\t.\t-\t.\tID=rna-Pt:301..420;gbkey=rRNA;product=5S ribosomal RNA\n"
        "Pt\tRefSeq\tgene\t501\t600\t.\t-\t.\tID=gene-L6;Name=L6;gene_biotype=other;locus_tag=L6\n"
        "Pt\tRefSeq\trRNA\t501\t600\t.\t-\t.\tID=rna-6;locus_tag=L7;product=4.5S ribosomal RNA\n"
        "Pt\tRefSeq\trRNA\t501\t600\t.\t+\t.\tID=rna-7;product=4.5S ribosomal RNA\n")
    genes = genes_from_gff3(read_gff3(gff)[0])
    by_name = {g.name: g for g in genes}
    assert [g.name for g in genes] == ["ATCG00010", "RRN16S", "LyesC2r005", "L6", "L7", "4.5S ribosomal RNA"]
    assert by_name["ATCG00010"].kind == "tRNA" and by_name["ATCG00010"].exons == [(4, 76)]
    assert by_name["RRN16S"].kind == "rRNA" and by_name["RRN16S"].symbol == "RRN16S"
    assert by_name["LyesC2r005"].kind == "rRNA" and by_name["LyesC2r005"].product == "5S ribosomal RNA"
    assert by_name["LyesC2r005"].label == "LyesC2r005 (5S ribosomal RNA)" and by_name["LyesC2r005"].symbol == ""
    assert by_name["L6"].kind == "other" and by_name["L7"].kind == "rRNA" and by_name["L7"].strand == -1
    assert by_name["4.5S ribosomal RNA"].strand == 1


def test_gff3_symbol_rules_one_by_one(tmp_path):
    gff = tmp_path / "symbols.gff3"
    gff.write_text(
        "##gff-version 3\n##sequence-region chr 1 1000\n"
        "chr\t.\tgene\t1\t30\t.\t+\t.\tID=g7;Name=L7;locus_tag=L7\n"  # The locus tag, whatever the ID
        "chr\t.\tCDS\t1\t30\t.\t+\t0\tID=c7;Parent=g7;product=protein seven\n"
        "chr\t.\tgene\t41\t60\t.\t+\t.\tID=g8;Name=AT8;gene_id=AT8\n"  # The gene_id (Ensembl)
        "chr\t.\tCDS\t41\t60\t.\t+\t0\tID=c8;Parent=g8;product=protein eight\n"
        "chr\t.\tCDS\t101\t130\t.\t+\t0\tID=c9;Name=YP_9;locus_tag=L9;product=protein nine\n"  # A CDS: never
        "chr\t.\tgene\t201\t230\t.\t+\t.\tID=g10;Name=psbA\n"  # A Name that is none of these: the symbol
        "chr\t.\tCDS\t201\t230\t.\t+\t0\tID=c10;Parent=g10;product=protein ten\n")
    genes = {g.name: g for g in genes_from_gff3(read_gff3(gff)[0])}
    assert genes["L7"].symbol == "" and genes["L7"].label == "L7 (protein seven)"
    assert genes["AT8"].symbol == "" and genes["AT8"].label == "AT8 (protein eight)"
    assert genes["L9"].symbol == "" and genes["L9"].label == "L9 (protein nine)"
    assert genes["psbA"].symbol == "psbA" and genes["psbA"].label == "psbA"


def test_diverged_stretches_are_aligned_in_reasonable_time():
    # Four 1.9 kb stretches of IRa with a mismatch every 20 bases (no 32-mer seed inside) and an indel each: the
    # stretches between seeds are aligned within a band (about 0.1 s; a full edit distance of two 1.9 kb
    # stretches takes about 2 s each, 14 s for the eight of the reviewer's case)
    import time

    from bacon.annotation import detect_inverted_repeat

    def edit(ira: str) -> str:
        s = list(ira)
        for block in range(4):
            start = 500 + block * 2500
            for p in range(start, start + 1900, 20):
                s[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[p]]
            s[start + 950] = ""
        return "".join(s)

    started = time.perf_counter()
    assert detect_inverted_repeat(make_plastome(ir=12000, edit=edit)) is None  # 380 differences in 12 kb
    assert time.perf_counter() - started < 3
