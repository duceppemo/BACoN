"""The bundled example (example/make_example.py): its annotated reference, planted SNPs and metadata agree with
what BACoN reads from them. The reads are not simulated here (that takes seconds); the reference is built by
importing the generator."""

import importlib.util
import random
import sys
from pathlib import Path

import pytest

from bacon.annotation import annotate_snps, genbank_fasta_records, load_annotation, translate
from bacon.metadata import choose_colour_column, read_metadata, unusable_reason
from bacon.seqio import Record, acgtn, write_fasta

EXAMPLE = Path(__file__).resolve().parent.parent / "example" / "make_example.py"


@pytest.fixture(scope="module")
def example():
    spec = importlib.util.spec_from_file_location("make_example", EXAMPLE)
    module = importlib.util.module_from_spec(spec)
    sys.modules["make_example"] = module  # dataclasses look the module up when the annotations are strings
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def reference(example, tmp_path_factory):
    out = tmp_path_factory.mktemp("example")
    ref = example.make_reference(random.Random(example.SEED))
    example.write_reference(out, ref)
    return out, ref


def test_reference_is_deterministic(example):
    a = example.make_reference(random.Random(example.SEED))
    b = example.make_reference(random.Random(example.SEED))
    assert a.seq == b.seq and len(a.seq) == example.GENOME_LENGTH
    assert [(s.position, s.alt) for s in a.snps] == [(s.position, s.alt) for s in b.snps]


def test_inverted_repeat_is_real(example, reference):
    _, ref = reference
    (b_start, b_end), (a_start, a_end) = example.REGIONS["IRb"], example.REGIONS["IRa"]
    assert ref.seq[a_start - 1:a_end] == example.revcomp(ref.seq[b_start - 1:b_end])


def test_genbank_gives_the_genes_and_regions(example, reference):
    out, ref = reference
    annotation = load_annotation(out / "reference.gb", [("organelle", example.GENOME_LENGTH)])
    assert annotation.warnings == [] and annotation.skipped == 0
    seq_ann = annotation.sequences["organelle"]
    assert seq_ann.circular
    assert [(r.name, r.start, r.end) for r in seq_ann.regions] == [(n, s, e) for n, (s, e) in example.REGIONS.items()]
    genes = {(g.name, g.kind, g.strand, tuple(g.exons or g.extent)) for g in seq_ann.genes}
    assert genes == {(f.name, f.kind, f.strand, tuple(f.parts)) for f in ref.features}
    kinds = [g.kind for g in seq_ann.genes]
    assert kinds.count("CDS") == 15 and kinds.count("pseudogene") == 1 and kinds.count("rRNA") == 2
    assert kinds.count("tRNA") == 5
    assert [g.name for g in seq_ann.genes if g.kind == "rRNA"] == ["rrn16-sim", "rrn16-sim"]  # Once per repeat
    split = next(g for g in seq_ann.genes if g.name == "orf04")
    assert len(split.cds[0].parts) == 2 and split.context(4851) == "intron"


def test_every_cds_translates_without_internal_stop(reference):
    out, ref = reference
    annotation = load_annotation(out / "reference.gb", [("organelle", len(ref.seq))])
    translations = {f.name: f.translation for f in ref.features if f.kind == "CDS"}
    for g in annotation.sequences["organelle"].genes:
        for cds in g.cds:
            coding = cds.coding_sequence(ref.seq)
            assert len(coding) % 3 == 0 and cds.transl_table == 11
            protein = "".join(translate(coding[i:i + 3], 11, start=i == 0) for i in range(0, len(coding), 3))
            assert protein[0] == "M" and protein[-1] == "*" and "*" not in protein[:-1], g.name
            assert protein[:-1] == translations[g.name]  # The /translation qualifier


def test_fasta_is_what_bacon_writes_from_the_genbank(reference, tmp_path):
    out, _ = reference
    records = [Record(r.header, acgtn(r.seq)) for r in genbank_fasta_records(out / "reference.gb")]
    assert [r.header for r in records] == ["organelle simulated circular reference"]
    write_fasta(tmp_path / "from_gb.fasta", records)
    assert (tmp_path / "from_gb.fasta").read_bytes() == (out / "reference.fasta").read_bytes()


def test_planted_snps_are_outside_the_repeats_and_apart(example, reference):
    out, ref = reference
    positions = [s.position for s in ref.snps]
    assert len(positions) == 20 and positions == sorted(positions)
    assert all(b - a >= example.MIN_SNP_SPACING for a, b in zip(positions, positions[1:]))
    assert positions[0] >= example.MIN_END_DISTANCE and positions[-1] <= len(ref.seq) - example.MIN_END_DISTANCE
    assert all(s.region in ("LSC", "SSC") for s in ref.snps)
    assert all(ref.seq[s.position - 1] == s.ref != s.alt for s in ref.snps)
    sets = ref.snp_sets()
    assert {k: len(v) for k, v in sets.items()} == {"b": 6, "g": 4, "d": 10}
    lines = (out / "planted_snps.tsv").read_text().splitlines()
    assert lines[0] == "sample\tpositions (1-based)" and lines[1] == "alpha\t"
    assert lines[3] == "gamma\t" + ",".join(str(p) for p, _ in sets["b"] + sets["g"])
    # The distances of the samples, from the sets they carry
    samples = {name: {c for k in keys for c in sets[k]} for name, keys in example.SAMPLES.items()}
    for (a, b), d in example.EXPECTED_DISTANCES.items():
        assert len(samples[a] ^ samples[b]) == d


def test_planted_effects_match_bacon(example, reference):
    out, ref = reference
    rows = [line.split("\t") for line in (out / "planted_effects.tsv").read_text().splitlines()]
    assert rows[0] == ["position", "ref", "alt", "samples", "region", "gene", "context", "codon", "amino_acid",
                       "effect"]
    annotation = load_annotation(out / "reference.gb", [("organelle", len(ref.seq))])
    info = annotate_snps(annotation, [("organelle", int(r[0]), r[1], r[2]) for r in rows[1:]], {"organelle": ref.seq})
    for r in rows[1:]:
        a = info[("organelle", int(r[0]))]
        found = [a.region, ", ".join(dict.fromkeys(g.name for g in a.genes)), a.context,
                 "; ".join(e.codons for e in a.effects), "; ".join(e.change for e in a.effects),
                 "; ".join(e.kind for e in a.effects)]
        assert found == r[4:10], r
    effects = [r[9] for r in rows[1:]]
    contexts = [r[6] for r in rows[1:]]
    assert {"synonymous", "missense", "nonsense"} <= set(effects)
    assert "intron" in contexts and "pseudogene" in contexts and contexts.count("tRNA") == 2
    assert sum(c.startswith("intergenic between") for c in contexts) == 4
    reverse = {f.name for f in ref.features if f.strand < 0}
    assert any(r[5] in reverse and r[9] for r in rows[1:])  # A coding change on the reverse strand
    assert any(r[5] == "orf04" and 5101 <= int(r[0]) <= 5699 for r in rows[1:])  # In the second exon


def test_metadata_has_the_samples_and_a_colour_column(example, reference):
    out, _ = reference
    assert (out / "metadata.tsv").read_text().startswith("# Simulated")
    metadata = read_metadata(out / "metadata.tsv")
    assert list(metadata.rows) == list(example.SAMPLES) == ["alpha", "beta", "gamma", "delta"]
    assert metadata.columns == ["group", "year", "origin", "note"]
    assert choose_colour_column(metadata, None) == ("group", None)
    assert metadata.value("delta", "group") == ""  # NA
    assert unusable_reason(metadata.values("note")) is not None  # Free text


CHECK = EXAMPLE.with_name("check_example.py")


@pytest.fixture(scope="module")
def checker(example):
    spec = importlib.util.spec_from_file_location("check_example", CHECK)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fake_run(example, reference, where):
    """An output folder of a perfect run: the truth's distances, a VCF of the planted SNPs and a report with the
    texts the check looks for."""
    out, ref = reference
    (where / "data").mkdir(parents=True)
    for name in ("reference.gb", "planted_effects.tsv"):
        (where / "data" / name).write_bytes((out / name).read_bytes())
    compared = where / "bacon" / "4_compared" / "ska"
    compared.mkdir(parents=True)
    (where / "bacon" / "reference.fasta").write_bytes((out / "reference.fasta").read_bytes())
    names = ["Reference", *example.SAMPLES]
    distances = {("Reference", s): d for (a, s), d in example.EXPECTED_DISTANCES.items() if a == "alpha"}
    distances.update({**example.EXPECTED_DISTANCES, ("Reference", "alpha"): 0})  # alpha is the reference
    distances.update({(b, a): d for (a, b), d in distances.items()})
    (compared / "snp_distances.tsv").write_text("snp-dists\t" + "\t".join(names) + "\n" + "".join(
        a + "".join(f"\t{0 if a == b else distances[(a, b)]}" for b in names) + "\n" for a in names))
    carried = {name: {c for k in keys for c in ref.snp_sets()[k]} for name, keys in example.SAMPLES.items()}
    (compared / "snps.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(example.SAMPLES)
        + "\n" + "".join(f"organelle\t{s.position}\t.\t{s.ref}\t{s.alt}\t.\t.\t.\tGT\t" + "\t".join(
            "1" if (s.position, s.alt) in carried[name] else "0" for name in example.SAMPLES) + "\n"
            for s in ref.snps))
    (where / "bacon" / "report.html").write_text(" ".join(text for text, _ in [
        ('<table class="snps', ""), ("in the LSC", ""), ("coloured by <b>group</b>", ""), ("site 2", "")]))
    return where


def test_check_passes_a_perfect_run(example, reference, checker, tmp_path, capsys):
    run = _fake_run(example, reference, tmp_path / "run")
    assert checker.check(run) == []
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 3 and all(line.startswith("OK: ") for line in lines)
    assert lines[1].startswith("OK: the 20 SNPs are annotated as planted (12 coding changes: ")


def test_check_fails_a_snp_on_another_sequence_or_a_wrong_effect(example, reference, checker, tmp_path, capsys):
    run = _fake_run(example, reference, tmp_path / "run")
    vcf = run / "bacon" / "4_compared" / "ska" / "snps.vcf"
    lines = vcf.read_text().splitlines()
    first = next(i for i, line in enumerate(lines) if not line.startswith("#"))
    lines[first] = "plasmid" + lines[first][len("organelle"):]  # The CHROM of one record
    vcf.write_text("\n".join(lines) + "\n")
    failures = checker.check(run)
    assert len(failures) == 1 and failures[0].startswith("VCF sites: ") and "('plasmid', " in failures[0]
    assert capsys.readouterr().out.count("OK: ") == 1  # The distances only
    # A wrong truth row: the annotation check names the SNP
    effects = run / "data" / "planted_effects.tsv"
    rows = effects.read_text().splitlines()
    cells = rows[2].split("\t")  # Not the first SNP, whose record is on the other sequence now
    cells[9] = "nonsense" if cells[9] != "nonsense" else "missense"
    rows[2] = "\t".join(cells)
    effects.write_text("\n".join(rows) + "\n")
    assert any(f.startswith(f"SNP {cells[0]}: ") for f in checker.check(run))


def test_check_main_exits_with_the_failures(example, reference, checker, tmp_path):
    run = _fake_run(example, reference, tmp_path / "run")
    (run / "bacon" / "report.html").write_text("nothing")
    with pytest.raises(SystemExit) as exc:
        checker.main(["check_example.py", str(run)])
    assert str(exc.value).startswith("FAILED: report.html lacks the SNP table")
    checker.main(["check_example.py", str(_fake_run(example, reference, tmp_path / "good"))])
