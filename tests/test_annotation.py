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
                   "c\tx\tgene\t70\t90\t.\t+\t.\tID=g3;gene_biotype=pseudogene;Name=g3\n")
    features, _, _ = read_gff3(gff)
    genes = genes_from_gff3(features)
    assert [(g.name, g.kind, g.strand) for g in genes] == [("psA", "pseudogene", 1), ("orphan", "CDS", -1),
                                                           ("g3", "pseudogene", 1)]
    assert genes[1].cds[0].codon_start == 3  # Phase 2 on the first part


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
    # Features beyond the end of the sequence are dropped with a warning
    ann = load_annotation(gb, [("ref.1", 50)])
    assert ann.genes == 1 and "5 feature(s) beyond the end" in ann.warnings[0] and ann.skipped == 5


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
                       ("inverted repeat region IRa", "a"), ("Inverted Repeat A region", "a"), ("IR_A", "a")]:
        match = _IR_COPY.search(text)
        assert match and match.group(1).lower() == copy, text
    for text in ("IRAK1 binding site", "inverted repeat", "spirAl", "LSC"):
        assert not _IR_COPY.search(text), text
