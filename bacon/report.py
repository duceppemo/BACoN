"""A self-contained HTML report of a BACoN output folder (standard library only).

Written at the end of every run as `report.html`; `python -m bacon.report OUTPUT` rebuilds it from the files of an
output folder (`run_info.json`, `summary.tsv`, the comparison's distances and tree).
"""

from __future__ import annotations

import html
import json
import math
import re
import sys
from pathlib import Path

from bacon.newick import parse
from bacon.seqio import split_extension

LOW_DEPTH = 20  # Same thresholds as the notes in summary.tsv
LENGTH_RANGE = (0.8, 1.2)
MAX_LABELLED_CELLS = 60  # Above this many genomes, heatmap cells show their value on hover only
MAX_HEATMAP_GENOMES = 150  # Above this, no inline heatmap (the page would be tens of MB): see the TSV

TABLE_COLUMNS = [
    ("Sample", "Sample"), ("Status", "Status"), ("Raw_reads", "Raw reads"), ("Baited_reads", "Baited reads"),
    ("Baited_pct", "Baited %"), ("Filtered_reads", "Filtered reads"), ("Filtered_N50", "Read N50"),
    ("Est_depth", "Depth"), ("Contigs", "Contigs"), ("Circular_contigs", "Circular"),
    ("Assembly_length", "Length"), ("Length_vs_reference", "× ref"), ("N_bases", "N bases"), ("Note", "Note"),
]
NUMERIC = {"Raw_reads", "Baited_reads", "Baited_pct", "Filtered_reads", "Filtered_N50", "Est_depth", "Contigs",
           "Circular_contigs", "Assembly_length", "Length_vs_reference", "N_bases"}


def _read_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line.strip()]


def _read_matrix(path: Path) -> tuple[list[str], dict[str, dict[str, int]]]:
    lines = path.read_text().splitlines()
    names = lines[0].split("\t")[1:]
    matrix = {}
    for line in lines[1:]:
        fields = line.split("\t")
        matrix[fields[0]] = dict(zip(names, map(int, fields[1:])))
    return names, matrix


def _float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _fmt(key: str, value: str) -> str:
    """Thousands separators for integers; values as written otherwise."""
    if key in NUMERIC and value.isdigit():
        return f"{int(value):,}"
    return value


def _flag(key: str, row: dict[str, str]) -> str:
    """CSS class for a summary cell that deserves attention."""
    value = row.get(key, "")
    if key == "Status" and value != "ok":
        return "bad"
    if key == "Est_depth" and (d := _float(value)) is not None and d < LOW_DEPTH:
        return "warn"
    if key == "Length_vs_reference" and (r := _float(value)) is not None and \
            not LENGTH_RANGE[0] <= r <= LENGTH_RANGE[1]:
        return "warn"
    if key == "N_bases" and value.isdigit() and int(value) > 0:
        return "info"
    return ""


def identical_groups(names: list[str], matrix: dict[str, dict[str, int]]) -> list[list[str]]:
    """Groups of two or more genomes with no SNP between any two of them.

    Distances ignore positions where either genome has N or a gap, so zero distances are not transitive
    (a = b and b = c does not imply a = c): a genome joins a group only if it is at distance 0 from every member.
    """
    groups: list[list[str]] = []
    for n in names:
        for g in groups:
            if all(matrix[n][m] == 0 for m in g):
                g.append(n)
                break
        else:
            groups.append([n])
    return sorted((g for g in groups if len(g) > 1), key=lambda g: (-len(g), g))


def distinct_count(names: list[str], matrix: dict[str, dict[str, int]]) -> int:
    groups = identical_groups(names, matrix)
    return len(groups) + sum(1 for n in names if not any(n in g for g in groups))


def _heat_color(t: float) -> tuple[str, str]:
    """Background and text colours for a value scaled to 0..1 (light yellow to dark blue)."""
    stops = [(255, 255, 217), (161, 218, 180), (65, 182, 196), (34, 94, 168), (8, 29, 88)]
    t = min(1.0, max(0.0, t))
    pos = t * (len(stops) - 1)
    i = min(int(pos), len(stops) - 2)
    f = pos - i
    rgb = [round(a + (b - a) * f) for a, b in zip(stops[i], stops[i + 1])]
    text = "#fff" if t > 0.55 else "#111"
    return f"rgb({rgb[0]},{rgb[1]},{rgb[2]})", text


def heatmap(names: list[str], matrix: dict[str, dict[str, int]]) -> str:
    """Distance matrix as an HTML table; colour on a log scale so that small distances stay visible."""
    top = max((matrix[a][b] for a in names for b in names), default=0)
    scale = math.log1p(top) or 1.0
    labelled = len(names) <= MAX_LABELLED_CELLS
    head = "".join(f'<th class="col"><span>{html.escape(n)}</span></th>' for n in names)
    rows = []
    for a in names:
        cells = []
        for b in names:
            d = matrix[a][b]
            bg, fg = _heat_color(math.log1p(d) / scale)
            title = f"{a} – {b}: {d} SNP{'s' if d != 1 else ''}"
            cells.append(f'<td style="background:{bg};color:{fg}" title="{html.escape(title)}">'
                         f'{d if labelled else ""}</td>')
        rows.append(f'<tr><th class="row">{html.escape(a)}</th>{"".join(cells)}</tr>')
    ticks = sorted({0, 1, 5, 10, 50, 100, 500, 1000, top} & set(range(top + 1))) if top else [0]
    legend = "".join(f'<span class="tick"><i style="background:{_heat_color(math.log1p(v) / scale)[0]}"></i>'
                     f"{v}</span>" for v in ticks)
    return (f'<div class="scroll"><table class="heat"><thead><tr><th></th>{head}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div><div class="legend">SNPs: {legend}</div>')


def _tool(info: dict, name: str) -> str:
    """' (version)' of a program from run_info.json: the version number only ('Filtlong v0.3.1' -> '0.3.1')."""
    text = (info.get("tools", {}).get(name) or {}).get("version") or ""
    match = re.search(r"\d+(?:\.\d+)+(?:-[\w.]+)?", text)
    return f" {html.escape(match.group(0))}" if match else ""


def _fasta_samples(info: dict) -> tuple[int, int]:
    """(samples with fasta reads, samples)."""
    samples = info.get("samples") or {}
    fasta = sum(1 for sample in samples.values()
                if any((split_extension(Path(f).name) or ("", ""))[1] == "fasta" for f in sample.get("files", [])))
    return fasta, len(samples)


def _num(value: object, default: float) -> float:
    return value if isinstance(value, (int, float)) and not isinstance(value, bool) else default


def methods_text(info: dict, rows: list[dict[str, str]] | None = None) -> str:
    """A methods paragraph built from the settings and program versions of the run."""
    s = info.get("settings") or {}
    ref = info.get("reference") or {}
    comparison = info.get("comparison") or {}
    parts = []
    ref_desc = f"{Path(ref.get('file') or s.get('reference') or 'the reference').name}"
    if ref.get("length"):
        ref_desc += f", {ref['length']:,} bp"
    if s.get("baiting") == "bbduk":
        parts.append(f"Reads sharing a {s.get('kmer', 31)}-mer (up to two mismatches) with the reference "
                     f"({ref_desc}) were extracted with BBDuk{_tool(info, 'bbduk.sh')}.")
    else:
        parts.append(f"Reads aligning to the reference ({ref_desc}) were extracted with "
                     f"minimap2{_tool(info, 'minimap2')} (-x map-ont).")
    size = "the given genome size" if s.get("genome_size") else "the reference length"
    min_length, depth = _num(s.get("min_read_length"), 500), _num(s.get("target_depth"), 100)
    fastq_text = (f"reads shorter than {min_length:g} bp were discarded and the best "
                  f"{_num(s.get('keep_percent'), 95):g}% were kept, up to {depth:g}x of {size}, with "
                  f"Filtlong{_tool(info, 'filtlong')}")
    fasta_text = (f"reads without qualities (fasta files) were filtered by BACoN: the longest reads of at "
                  f"least {min_length:g} bp, up to {depth:g}x of {size}")
    fasta, total = _fasta_samples(info)
    if fasta == 0:
        parts.append(fastq_text[0].upper() + fastq_text[1:] + ".")
    elif fasta == total:
        parts.append(fasta_text[0].upper() + fasta_text[1:] + ".")
    else:
        parts.append(f"For samples with base qualities, {fastq_text}; {fasta_text}.")
    assembler = s.get("assembler")
    if assembler == "samtools":
        gaps = ("; runs of N at the ends of a sequence were replaced by the reference"
                if s.get("template_gaps") == "reference" else "")
        parts.append("Each sample was assembled by reference-guided consensus: reads were aligned to the "
                     "reference with minimap2 and the consensus called with samtools consensus"
                     f"{_tool(info, 'samtools')} (-X r10.4_sup, minimum depth 3; positions with less support, or "
                     f"where the reads disagree, are N{gaps}).")
    elif assembler == "flye":
        overlap = f", minimum overlap {s['min_size']} bp" if s.get("min_size") else ""
        parts.append(f"Each sample was assembled de novo with Flye{_tool(info, 'flye')} (--{s.get('read_type')}, "
                     f"genome size {'as given' if s.get('genome_size') else 'the reference length'}, "
                     f"{s.get('flye_iterations')} polishing iterations{overlap}).")
    elif assembler == "myloasm":
        parts.append(f"Each sample was assembled de novo with myloasm{_tool(info, 'myloasm')}.")
    added = s.get("add_genomes") or []
    if added:
        parts.append(f"{len(added)} finished genome{'s were' if len(added) > 1 else ' was'} added to the "
                     "comparison.")
    method = s.get("snp_method")
    compared = comparison.get("distances") and not comparison.get("failed") and not comparison.get("skipped")
    circular = sum(int(r["Circular_contigs"]) for r in rows or [] if r.get("Circular_contigs", "").isdigit())
    if method == "ska" and compared:
        freq = _num(s.get("ska_min_freq"), 1.0)
        scope = ("present in all genomes (core SNPs)" if freq >= 1
                 else f"present in at least {freq * 100:g}% of the genomes")
        wrap = (", circular contigs (and the reference, when most assemblies were circular) being extended by 30 "
                "bases so that variants near their ends are kept" if circular else "")
        parts.append(f"SNPs were identified with SKA2{_tool(info, 'ska')} from split 31-mers {scope}{wrap}.")
    elif method == "parsnp" and compared:
        parts.append(f"The core genome was aligned to the reference with Parsnp{_tool(info, 'parsnp')} (-c); SNPs "
                     f"were extracted with HarvestTools{_tool(info, 'harvesttools')}.")
    if compared:
        parts.append("Pairwise SNP distances count the positions where both genomes have a nucleotide and they "
                     "differ.")
    if compared and comparison.get("vcf"):
        how = (f"by mapping the split k-mers to the reference with ska map{_tool(info, 'ska')}"
               if method == "ska" else f"from the Parsnp alignment with HarvestTools{_tool(info, 'harvesttools')}")
        parts.append(f"The SNPs of each genome relative to the reference were written to a VCF file {how}.")
    if comparison.get("tree"):
        on = "the core-genome alignment" if method == "parsnp" else "the SNP alignment"
        if s.get("tree") == "iqtree":
            parts.append(f"A tree was built on {on} with IQ-TREE{_tool(info, 'iqtree')} (ModelFinder, 1000 "
                         "ultrafast bootstraps; the bootstrap consensus tree is shown) and rooted at its "
                         "midpoint.")
        else:
            parts.append(f"A tree was built on {on} with FastTree{_tool(info, 'FastTree')} (GTR, SH-like supports "
                         "from 100 resamples) and rooted at its midpoint.")
    parts.append(f"The analysis was run with BACoN {html.escape(str(info.get('bacon_version', '')))} "
                 "(https://github.com/duceppemo/BACoN).")
    return " ".join(parts)


CSS = """
:root {
  --bg: #fff; --fg: #1d2433; --muted: #5d6678; --line: #e3e6ec; --card: #f6f7f9; --accent: #b63a2f;
  --bad: #fde2e1; --bad-fg: #9b1c1c; --warn: #fdf0d5; --warn-fg: #7a4b00; --info: #e3eefc; --info-fg: #1e4f91;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #15181e; --fg: #e6e8ec; --muted: #9aa3b2; --line: #2c313b; --card: #1d2129; --accent: #f0806f;
    --bad: #4a1f1f; --bad-fg: #ffb4ab; --warn: #45371a; --warn-fg: #f5cf82; --info: #1c3050; --info-fg: #a9c8f5;
  }
}
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--bg); color: var(--fg);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
}
main { max-width: 1200px; margin: 0 auto; padding: 24px 16px 64px; }
h1 { font-size: 26px; margin: 0 0 4px; }
h1 b { color: var(--accent); }
h2 { font-size: 19px; margin: 36px 0 10px; padding-top: 8px; border-top: 1px solid var(--line); }
.sub { color: var(--muted); margin: 0 0 20px; }
.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 10px; }
.card { background: var(--card); border-radius: 8px; padding: 10px 14px; }
.card b { display: block; font-size: 22px; }
.card span { color: var(--muted); font-size: 13px; }
.scroll { overflow-x: auto; max-width: 100%; }
table { border-collapse: collapse; font-size: 13px; }
th, td { padding: 4px 8px; text-align: left; white-space: nowrap; }
table.samples th {
  position: sticky; top: 0; background: var(--bg); cursor: pointer; user-select: none;
  border-bottom: 2px solid var(--line);
}
table.samples th:hover { color: var(--accent); }
table.samples td { border-bottom: 1px solid var(--line); }
table.samples td.num { text-align: right; font-variant-numeric: tabular-nums; }
table.samples td.note { white-space: normal; min-width: 220px; }
.bad { background: var(--bad); color: var(--bad-fg); }
.warn { background: var(--warn); color: var(--warn-fg); }
.info { background: var(--info); color: var(--info-fg); }
.key span { display: inline-block; padding: 1px 8px; border-radius: 4px; margin-right: 8px; font-size: 13px; }
.tree { background: #fff; border-radius: 8px; padding: 8px; overflow-x: auto; }
.tree svg { max-width: 100%; height: auto; }
table.heat td {
  width: 26px; min-width: 26px; text-align: center; padding: 2px; font-size: 11px;
  font-variant-numeric: tabular-nums;
}
table.heat th.row { text-align: right; font-weight: normal; font-size: 12px; }
table.heat th.col { height: 120px; vertical-align: bottom; padding: 0; }
table.heat th.col span {
  display: inline-block; writing-mode: vertical-rl; transform: rotate(180deg); font-weight: normal;
  font-size: 12px; padding: 4px 0;
}
.legend { margin-top: 8px; color: var(--muted); font-size: 13px; }
.legend .tick { margin-right: 10px; }
.legend i {
  display: inline-block; width: 14px; height: 14px; vertical-align: -2px; margin-right: 4px;
  border: 1px solid var(--line);
}
ul.groups li { margin-bottom: 4px; }
.methods { background: var(--card); border-radius: 8px; padding: 12px 16px; }
dl { display: grid; grid-template-columns: max-content 1fr; gap: 4px 16px; font-size: 13px; }
dt { color: var(--muted); }
dd { margin: 0; word-break: break-all; }
code { font-size: 12px; }
table.heat td, .bad, .warn, .info, .legend i { print-color-adjust: exact; -webkit-print-color-adjust: exact; }
@media print {
  .tree { border: 1px solid #ccc; }
  table.samples th { position: static; }
}
"""

SORT_JS = """
document.querySelectorAll('table.samples th').forEach(function (th, i) {
  th.addEventListener('click', function () {
    var table = th.closest('table'), body = table.tBodies[0], rows = Array.from(body.rows);
    var ascending = th.dataset.asc !== '1', numeric = th.classList.contains('num');
    table.querySelectorAll('th').forEach(function (h) { h.dataset.asc = ''; });
    th.dataset.asc = ascending ? '1' : '0';
    rows.sort(function (x, y) {
      var a = x.cells[i].dataset.v, b = y.cells[i].dataset.v;
      var r = numeric ? parseFloat(a) - parseFloat(b) : a.localeCompare(b, undefined, {numeric: true});
      return ascending ? r : -r;
    });
    rows.forEach(function (r) { body.appendChild(r); });
  });
});
"""


def build_report(output: Path) -> str:
    info = json.loads((output / "run_info.json").read_text())
    rows = _read_tsv(output / "summary.tsv") if (output / "summary.tsv").exists() else []
    comparison = info.get("comparison") or {}
    ref = info.get("reference") or {}
    ok = sum(1 for r in rows if r.get("Status") == "ok")
    flagged = sum(1 for r in rows if r.get("Status") == "ok" and r.get("Note"))
    esc = html.escape

    out = ["<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\">",
           '<meta name="viewport" content="width=device-width,initial-scale=1">',
           f"<title>BACoN report</title><style>{CSS}</style></head><body><main>",
           "<h1><b>BACoN</b> report</h1>",
           f'<p class="sub">{esc(Path(str(output)).name)} · {esc(str(info.get("started", "")))} · '
           f'BACoN {esc(str(info.get("bacon_version", "")))} · last run {info.get("duration_s", "?")} s</p>']

    # Overview
    cards = [(f"{ok}/{len(rows)}", "samples assembled"), (str(flagged), "assembled with a note"),
             (f"{ref.get('length', 0):,} bp" if ref.get("length") else "?", "reference"),
             (str(comparison.get("core_snps", "–")), "SNP sites")]
    out.append('<div class="cards">' + "".join(f'<div class="card"><b>{esc(v)}</b><span>{esc(k)}</span></div>'
                                               for v, k in cards) + "</div>")
    if comparison.get("failed") and isinstance(comparison["failed"], str):
        out.append(f'<p class="bad">The comparison failed: {esc(comparison["failed"])}</p>')
    elif comparison.get("skipped"):
        out.append(f"<p>The comparison was skipped: {esc(str(comparison['skipped']))}.</p>")

    # Samples
    out.append("<h2>Samples</h2>")
    out.append('<p class="key"><span class="bad">failed</span>'
               f'<span class="warn">depth below {LOW_DEPTH}x or length outside '
               f'{LENGTH_RANGE[0]}–{LENGTH_RANGE[1]}x the reference</span>'
               '<span class="info">N bases</span>Click a column to sort.</p>')
    out.append('<div class="scroll"><table class="samples"><thead><tr>'
               + "".join(f'<th class="{"num" if key in NUMERIC else ""}">{esc(label)}</th>'
                         for key, label in TABLE_COLUMNS) + "</tr></thead><tbody>")
    for r in rows:
        cells = []
        for key, _ in TABLE_COLUMNS:
            value = r.get(key, "")
            classes = " ".join(c for c in (_flag(key, r), "num" if key in NUMERIC else "",
                                           "note" if key == "Note" else "") if c)
            sort_value = value if value not in ("NA", "") else ("-1" if key in NUMERIC else "")
            cells.append(f'<td class="{classes}" data-v="{esc(sort_value)}">{esc(_fmt(key, value))}</td>')
        out.append("<tr>" + "".join(cells) + "</tr>")
    out.append("</tbody></table></div>")

    # Comparison
    notes: list[str] = []

    def located(key: str) -> Path | None:
        """A file of the comparison: relative paths start from the output folder; an absolute path that no
        longer exists (the folder was moved) is looked up under 4_compared/ of this folder."""
        if not comparison.get(key):
            return None
        path = output / comparison[key]
        if not path.exists():
            parts = Path(comparison[key]).parts
            if "4_compared" in parts:
                moved = output.joinpath(*parts[parts.index("4_compared"):])
                if moved.exists():
                    return moved
            notes.append(f"{comparison[key]} was not found")
        return path

    distances = located("distances")
    if distances is not None and distances.exists():
        names, matrix = _read_matrix(distances)
        tree_file = located("tree")
        order, in_tree_order = names, False
        if tree_file is not None and tree_file.exists():
            leaves = [leaf.name for leaf in parse(tree_file.read_text()).leaves()]
            if sorted(leaves) == sorted(names):
                order, in_tree_order = leaves, True
        method = comparison.get("method", "")
        method_name = {"ska": "SKA2", "parsnp": "Parsnp"}.get(method, method)
        tree_tool = {"fasttree": "FastTree", "iqtree": "IQ-TREE"}.get(comparison.get("tree_method", ""), "")
        svg = located("tree_svg")
        out.append("<h2>Tree</h2>")
        if svg is not None and svg.exists():
            out.append(f'<p class="sub">{esc(method_name)} SNPs; {esc(tree_tool)}, midpoint-rooted; internal '
                       "labels are supports.</p>")
            out.append(f'<div class="tree">{svg.read_text()}</div>')
        elif comparison.get("core_snps") == 0:
            why = ("no SNP site is shared by all the genomes; see --ska-min-freq"
                   if method == "ska" else "the alignment has no SNP site")
            out.append(f"<p>No tree: {why}.</p>")
        else:
            out.append("<p>No tree.</p>")
        out.append("<h2>SNP distances</h2>")
        order_text = "in tree order" if in_tree_order else "in the order of the distance table"
        out.append(f'<p class="sub">Pairwise SNP distances, {order_text}. Hover a cell for its pair.</p>')
        if len(order) <= MAX_HEATMAP_GENOMES:
            out.append(heatmap(order, matrix))
        else:
            shown = distances.relative_to(output) if distances.is_relative_to(output) else distances
            out.append(f"<p>{len(order)} genomes: too many for a heatmap in this page; the distances are in "
                       f"<code>{esc(str(shown))}</code>.</p>")
        if comparison.get("core_snps") == 0:
            out.append("<p>No SNP site was compared: the distances say nothing about identity.</p>")
        else:
            groups = identical_groups(order, matrix)
            if groups:
                out.append("<h3>Identical genomes</h3><p class=\"sub\">No SNP between any two genomes of a group "
                           "(positions with N or a gap are not compared).</p><ul class=\"groups\">" + "".join(
                               f"<li><b>{len(g)}</b>: {esc(', '.join(g))}</li>" for g in groups) + "</ul>")
            out.append(f"<p>{len(order)} genomes, {distinct_count(order, matrix)} distinct at the SNP sites "
                       "compared.</p>")
    for note in notes:
        out.append(f'<p class="warn">{esc(note)}.</p>')

    # Methods and provenance
    out.append(f'<h2>Methods</h2><p class="methods">{methods_text(info, rows)}</p>')
    out.append("<h2>Run</h2><dl>")
    prov = [("Command", " ".join(info.get("command_line") or [])),
            ("Reference", f"{ref.get('file', '')} ({ref.get('sequences', '?')} sequence(s), "
                          f"MD5 {ref.get('md5', '?')})"),
            ("Output", str(output)), ("Python", f"{info.get('python', '')} on {info.get('platform', '')}")]
    prov += [(name, f"{(t or {}).get('version') or '?'} — {(t or {}).get('path', '')}")
             for name, t in (info.get("tools") or {}).items()]
    out.append("".join(f"<dt>{esc(k)}</dt><dd><code>{esc(v)}</code></dd>" for k, v in prov))
    out.append("</dl></main>")
    out.append(f"<script>{SORT_JS}</script></body></html>")
    return "\n".join(out) + "\n"


def write_report(output: Path) -> Path:
    path = output / "report.html"
    tmp = path.with_suffix(".html.tmp")
    tmp.write_text(build_report(output), encoding="utf-8")
    tmp.replace(path)
    return path


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python -m bacon.report OUTPUT_FOLDER")
    print(write_report(Path(sys.argv[1])))
