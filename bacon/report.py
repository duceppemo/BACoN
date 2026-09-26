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
    """Groups of two or more genomes at distance 0 from one another (connected through zero distances)."""
    parent = {n: n for n in names}

    def find(n: str) -> str:
        while parent[n] != n:
            parent[n] = parent[parent[n]]
            n = parent[n]
        return n

    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if matrix[a][b] == 0:
                parent[find(a)] = find(b)
    groups: dict[str, list[str]] = {}
    for n in names:
        groups.setdefault(find(n), []).append(n)
    return sorted((g for g in groups.values() if len(g) > 1), key=lambda g: (-len(g), g))


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


def _fasta_input(info: dict) -> bool:
    files = [f for sample in info.get("samples", {}).values() for f in sample.get("files", [])]
    return any((split_extension(Path(f).name) or ("", ""))[1] == "fasta" for f in files)


def methods_text(info: dict) -> str:
    """A methods paragraph built from the settings and program versions of the run."""
    s = info.get("settings", {})
    ref = info.get("reference") or {}
    parts = []
    ref_desc = f"{Path(ref.get('file') or s.get('reference', 'the reference')).name}"
    if ref.get("length"):
        ref_desc += f", {ref['length']:,} bp"
    if s.get("baiting") == "bbduk":
        parts.append(f"Reads sharing a {s.get('kmer')}-mer (up to two mismatches) with the reference ({ref_desc}) "
                     f"were extracted with BBDuk{_tool(info, 'bbduk.sh')}.")
    else:
        parts.append(f"Reads aligning to the reference ({ref_desc}) were extracted with "
                     f"minimap2{_tool(info, 'minimap2')} (-x map-ont).")
    parts.append(f"Reads shorter than {s.get('min_read_length')} bp were discarded and the best "
                 f"{s.get('keep_percent'):g}% were kept, up to {s.get('target_depth')}x of the "
                 f"{'given genome size' if s.get('genome_size') else 'reference length'}, with "
                 f"Filtlong{_tool(info, 'filtlong')}"
                 + (" (reads from fasta files, without qualities: the longest ones, selected by BACoN)."
                    if _fasta_input(info) else "."))
    assembler = s.get("assembler")
    if assembler == "samtools":
        parts.append("Each sample was assembled by reference-guided consensus: reads were aligned to the "
                     "reference with minimap2 and the consensus called with samtools consensus"
                     f"{_tool(info, 'samtools')} "
                     "(-X r10.4_sup, minimum depth 3; positions with less support are N).")
    elif assembler == "flye":
        parts.append(f"Each sample was assembled de novo with Flye{_tool(info, 'flye')} "
                     f"(--{s.get('read_type')}, {s.get('flye_iterations')} polishing iterations).")
    elif assembler == "myloasm":
        parts.append(f"Each sample was assembled de novo with myloasm{_tool(info, 'myloasm')}.")
    method = s.get("snp_method")
    comparison = info.get("comparison") or {}
    if method == "ska" and not comparison.get("skipped"):
        freq = s.get("ska_min_freq", 1)
        scope = "present in all genomes (core SNPs)" if freq == 1 else f"present in at least {freq:.0%} of genomes"
        wrap = (", circular contigs being extended by 30 bases so that variants near their ends are kept"
                if assembler in ("flye", "myloasm") else "")
        parts.append(f"SNPs were identified with SKA2{_tool(info, 'ska')} from split 31-mers {scope}{wrap}.")
    elif method == "parsnp" and not comparison.get("skipped"):
        parts.append(f"Core-genome SNPs were identified with Parsnp{_tool(info, 'parsnp')} (-c) and "
                     f"HarvestTools{_tool(info, 'harvesttools')}.")
    if comparison.get("tree"):
        if s.get("tree") == "iqtree":
            parts.append(f"A tree was built with IQ-TREE{_tool(info, 'iqtree')} (ModelFinder, 1000 ultrafast "
                         "bootstraps) and rooted at its midpoint.")
        else:
            parts.append(f"A tree was built with FastTree{_tool(info, 'FastTree')} (GTR, SH-like supports from "
                         "100 resamples) and rooted at its midpoint.")
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
@media print {
  .tree { border: 1px solid #ccc; }
  table.samples th { position: static; }
}
"""

SORT_JS = """
document.querySelectorAll('table.samples th').forEach(function (th, i) {
  th.addEventListener('click', function () {
    var table = th.closest('table'), body = table.tBodies[0], rows = Array.from(body.rows);
    var ascending = th.dataset.asc !== '1';
    table.querySelectorAll('th').forEach(function (h) { h.dataset.asc = ''; });
    th.dataset.asc = ascending ? '1' : '0';
    rows.sort(function (x, y) {
      var a = x.cells[i].dataset.v, b = y.cells[i].dataset.v, na = parseFloat(a), nb = parseFloat(b);
      var r = (!isNaN(na) && !isNaN(nb)) ? na - nb : a.localeCompare(b);
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
               + "".join(f"<th>{esc(label)}</th>" for _, label in TABLE_COLUMNS) + "</tr></thead><tbody>")
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
    distances = Path(comparison["distances"]) if comparison.get("distances") else None
    if distances is not None and not distances.is_absolute():
        distances = output / distances
    if distances is not None and distances.exists():
        names, matrix = _read_matrix(distances)
        tree_file = Path(comparison["tree"]) if comparison.get("tree") else None
        order = names
        if tree_file is not None and tree_file.exists():
            leaves = [leaf.name for leaf in parse(tree_file.read_text()).leaves()]
            if sorted(leaves) == sorted(names):
                order = leaves
        method = comparison.get("method", "")
        method = {"ska": "SKA2", "parsnp": "Parsnp"}.get(method, method)
        out.append(f"<h2>Tree</h2><p class=\"sub\">{esc(method)} SNPs; "
                   f"{'FastTree' if comparison.get('tree_method') == 'fasttree' else 'IQ-TREE'}, "
                   "midpoint-rooted; internal labels are supports.</p>")
        svg = Path(comparison["tree_svg"]) if comparison.get("tree_svg") else None
        if svg is not None and svg.exists():
            out.append(f'<div class="tree">{svg.read_text()}</div>')
        else:
            out.append("<p>No tree: no SNP site is shared by all the genomes.</p>")
        out.append("<h2>SNP distances</h2>"
                   '<p class="sub">Pairwise SNP distances, in tree order. Hover a cell for its pair.</p>')
        out.append(heatmap(order, matrix))
        groups = identical_groups(order, matrix)
        if groups:
            out.append("<h3>Identical genomes</h3><ul class=\"groups\">" + "".join(
                f"<li><b>{len(g)}</b>: {esc(', '.join(g))}</li>" for g in groups) + "</ul>")
        distinct = len(identical_groups(order, matrix)) + sum(
            1 for n in order if not any(n in g for g in groups))
        out.append(f"<p>{len(order)} genomes, {distinct} distinct at the SNP sites compared.</p>")

    # Methods and provenance
    out.append(f'<h2>Methods</h2><p class="methods">{methods_text(info)}</p>')
    out.append("<h2>Run</h2><dl>")
    prov = [("Command", " ".join(info.get("command_line", []))),
            ("Reference", f"{ref.get('file', '')} ({ref.get('sequences', '?')} sequence(s), "
                          f"MD5 {ref.get('md5', '?')})"),
            ("Output", str(output)), ("Python", f"{info.get('python', '')} on {info.get('platform', '')}")]
    prov += [(name, f"{t.get('version') or '?'} — {t.get('path', '')}")
             for name, t in info.get("tools", {}).items()]
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
