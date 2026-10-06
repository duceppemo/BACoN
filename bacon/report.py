"""A self-contained HTML report of a BACoN output folder (standard library only).

Written at the end of every run as `report.html`; `python -m bacon.report OUTPUT` rebuilds it from the files of an
output folder (`run_info.json`, `summary.tsv`, `reference.fasta`, the assemblies and the comparison's distances,
tree and VCF). Every figure is inline SVG drawn here; a missing or unreadable file only removes its figure.
"""

from __future__ import annotations

import html
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

from bacon.newick import Node, ladderize, parse
from bacon.seqio import read_records, split_extension

LOW_DEPTH = 20  # Same thresholds as the notes in summary.tsv
LENGTH_RANGE = (0.8, 1.2)
MAX_LABELLED_CELLS = 60  # Above this many genomes, heatmap cells show their value on hover only
MAX_HEATMAP_GENOMES = 150  # Above this, no inline heatmap (the page would be tens of MB): see the TSV
MAX_MAP_SEQUENCES = 8  # Reference sequences drawn in the genome map (the longest ones)
GROUP_COLOURS = 8  # Groups of identical genomes with a colour of their own (--s1 to --s8); the others share grey
N_LABELS = 5  # Bars labelled in the N-bases chart (the largest values)
FIG_WIDTH = 1120

TABLE_COLUMNS = [
    ("Sample", "Sample"), ("Status", "Status"), ("Raw_reads", "Raw reads"), ("Baited_reads", "Baited reads"),
    ("Baited_pct", "Baited %"), ("Filtered_reads", "Filtered reads"), ("Filtered_N50", "Read N50"),
    ("Est_depth", "Depth"), ("Contigs", "Contigs"), ("Circular_contigs", "Circular"),
    ("Assembly_length", "Length"), ("Length_vs_reference", "× ref"), ("N_bases", "N bases"), ("Note", "Note"),
]
NUMERIC = {"Raw_reads", "Baited_reads", "Baited_pct", "Filtered_reads", "Filtered_N50", "Est_depth", "Contigs",
           "Circular_contigs", "Assembly_length", "Length_vs_reference", "N_bases"}

esc = html.escape


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


def group_slots(groups: list[list[str]]) -> dict[str, int]:
    """Genome -> index of its group of identical genomes (groups as returned by identical_groups)."""
    return {name: i for i, g in enumerate(groups) for name in g}


def slot_colour(slot: int) -> str:
    return f"var(--s{slot + 1})" if slot < GROUP_COLOURS else "var(--muted)"


# ---------------------------------------------------------------------------------------------------------------
# SVG helpers
# ---------------------------------------------------------------------------------------------------------------

def _svg(width: float, height: float, title: str) -> str:
    return (f'<svg class="fig" viewBox="0 0 {width:.0f} {height:.0f}" width="100%" '
            f'style="max-width:{width:.0f}px" role="img"><title>{esc(title)}</title>')


def _text_px(text: str, size: float = 11) -> float:
    """Rough width of a string in the system sans: enough to reserve room for labels."""
    return 0.58 * size * len(text)


def _nice_ceil(x: float) -> float:
    """The smallest of 1, 2, 5 x 10^k at or above x."""
    if x <= 0:
        return 1.0
    exponent = 10 ** math.floor(math.log10(x))
    for step in (1, 2, 5, 10):
        if step * exponent >= x:
            return step * exponent
    return 10 * exponent  # pragma: no cover


def _nice_floor(x: float) -> int:
    """The largest of 1, 2, 3, 5 x 10^k at or below x (x >= 1)."""
    exponent = 10 ** math.floor(math.log10(max(x, 1.0)))
    for step in (5, 3, 2, 1):
        if step * exponent <= x:
            return int(step * exponent)
    return 1  # pragma: no cover


def _ticks(top: float, count: int = 5, minimum_step: float = 0.0) -> list[float]:
    step = max(_nice_ceil(top / count) if top > 0 else 1.0, minimum_step)
    return [i * step for i in range(int(top // step) + 1)]


def _g(x: float) -> str:
    return f"{x:g}"


# ---------------------------------------------------------------------------------------------------------------
# Figure 1 and 2: per-sample bars
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class Bar:
    name: str
    value: float | None  # None: no value (NA, or a failed sample)
    hover: str
    label: str = ""  # Text at the end of the bar (flagged samples only)
    failed: str = ""  # The status of a failed sample ("failed (bait)"): listed without a bar


def bar_chart(bars: list[Bar], title: str, *, unit: str = "", threshold: float | None = None,
              threshold_text: str = "", width: int = 560, integers: bool = False) -> str:
    """Horizontal bars sorted by value, largest first; a threshold line; failed samples listed at the end."""
    bars = (sorted((b for b in bars if b.value is not None and not b.failed), key=lambda b: -b.value)
            + [b for b in bars if b.value is None and not b.failed] + [b for b in bars if b.failed])
    row, size = 16, 11
    left = min(220, 16 + _text_px(max((b.name for b in bars), key=len, default=""), size))
    plot = width - left - 80
    top = 26
    values = [b.value for b in bars if b.value is not None]
    x_max = max(values + ([threshold] if threshold else [])) if values or threshold else 0
    x_max = x_max * 1.05 if x_max > 0 else 1.0
    sx = plot / x_max
    height = top + len(bars) * row + (34 if threshold else 16)
    parts = [_svg(width, height, title)]
    for t in _ticks(x_max, minimum_step=1.0 if integers else 0.0):
        x = left + t * sx
        parts.append(f'<line x1="{x:.1f}" y1="{top - 4}" x2="{x:.1f}" y2="{top + len(bars) * row}" '
                     'class="grid"/>')
        parts.append(f'<text x="{x:.1f}" y="{top - 8}" text-anchor="middle" class="t-tiny t-muted">'
                     f'{_g(t)}{unit}</text>')
    if threshold is not None:
        x = left + threshold * sx
        parts.append(f'<line x1="{x:.1f}" y1="{top - 4}" x2="{x:.1f}" y2="{top + len(bars) * row}" '
                     'class="thresh"/>')
        parts.append(f'<text x="{x + 4:.1f}" y="{top + len(bars) * row + 14}" class="t-tiny t-ink2">'
                     f'{esc(threshold_text)}</text>')
    for i, b in enumerate(bars):
        y = top + i * row
        parts.append(f'<text x="{left - 6}" y="{y + row - 4}" text-anchor="end" class="t-small t-ink2">'
                     f'{esc(b.name)}</text>')
        if b.failed:
            parts.append(f'<text x="{left}" y="{y + row - 4}" class="t-tiny t-muted"><tspan fill="var(--crit)">'
                         f'&#9632;</tspan> {esc(b.failed)}<title>{esc(b.hover)}</title></text>')
        elif b.value is None:
            parts.append(f'<text x="{left}" y="{y + row - 4}" class="t-tiny t-muted">NA<title>{esc(b.hover)}'
                         '</title></text>')
        else:
            w = max(0.0, b.value * sx)
            cls = "bar" if w >= 1 else "bar0"
            parts.append(f'<rect x="{left}" y="{y + 3}" width="{max(w, 1.0):.1f}" height="{row - 6}" rx="2" '
                         f'class="{cls}"><title>{esc(b.hover)}</title></rect>')
            if b.label:
                parts.append(f'<text x="{left + w + 5:.1f}" y="{y + row - 4}" class="t-tiny t-ink">'
                             f'{esc(b.label)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def depth_bars(rows: list[dict[str, str]]) -> list[Bar]:
    bars = []
    for r in rows:
        name, status = r.get("Sample", ""), r.get("Status", "")
        depth = _float(r.get("Est_depth", ""))
        hover = f"{name}: {f'{depth:g}x' if depth is not None else 'no'} depth after filtering"
        if r.get("Filtered_reads", "").isdigit():
            n50 = _fmt("Filtered_N50", r.get("Filtered_N50", "?"))
            hover += f"; {int(r['Filtered_reads']):,} reads, N50 {n50}"
        if r.get("Note"):
            hover += f". {r['Note']}"
        flagged = _flag("Est_depth", r) == "warn"
        bars.append(Bar(name, depth, hover, label=f"{depth:g}x" if flagged and depth is not None else "",
                        failed="" if status == "ok" else status or "failed"))
    return bars


def n_bars(rows: list[dict[str, str]]) -> list[Bar]:
    bars = []
    values = sorted((int(r["N_bases"]) for r in rows if r.get("N_bases", "").isdigit()), reverse=True)
    label_from = values[N_LABELS - 1] if len(values) >= N_LABELS else (values[-1] if values else 0)
    for r in rows:
        name, status = r.get("Sample", ""), r.get("Status", "")
        n = int(r["N_bases"]) if r.get("N_bases", "").isdigit() else None
        hover = f"{name}: {n:,} N bases" if n is not None else f"{name}: no N count"
        if n is not None and r.get("Assembly_length", "").isdigit() and int(r["Assembly_length"]):
            hover += f" of {int(r['Assembly_length']):,} bp ({100 * n / int(r['Assembly_length']):.2f}%)"
        if r.get("Note"):
            hover += f". {r['Note']}"
        label = f"{n:,}" if n is not None and n > 0 and n >= label_from else ""
        bars.append(Bar(name, n, hover, label=label, failed="" if status == "ok" else status or "failed"))
    return bars


# ---------------------------------------------------------------------------------------------------------------
# Figure: the tree
# ---------------------------------------------------------------------------------------------------------------

def tree_svg(root: Node, slots: dict[str, int], ref_name: str = "", snp_sites: int | None = None) -> str:
    """A rectangular phylogram: leaves with a square in the colour of their group of identical genomes, the
    Reference in bold, supports on the internal branches, a scale bar in substitutions per site."""
    leaves = root.leaves()
    row, size = (18, 12) if len(leaves) <= 80 else (14, 10.5)
    left, top, plot = 16.0, 14.0, 560.0
    depth: dict[int, float] = {}

    def assign_x(n: Node, d: float) -> None:
        depth[id(n)] = d
        for c in n.children:
            assign_x(c, d + max(c.length, 0.0))

    assign_x(root, 0.0)
    max_depth = max(depth.values())
    sx = plot / max_depth if max_depth > 0 else 0.0
    y: dict[int, float] = {id(leaf): top + i * row + row / 2 for i, leaf in enumerate(leaves)}

    def assign_y(n: Node) -> float:
        if n.is_leaf():
            return y[id(n)]
        ys = [assign_y(c) for c in n.children]
        y[id(n)] = (min(ys) + max(ys)) / 2
        return y[id(n)]

    assign_y(root)
    longest = max((len(leaf.name) + (len(ref_name) + 1 if leaf.name == "Reference" else 0) for leaf in leaves),
                  default=1)
    width = left + plot + 30 + _text_px("x" * longest, size) + (13 if slots else 0)
    height = top + len(leaves) * row + (40 if max_depth > 0 else 12)
    parts = [_svg(width, height, "Tree")]
    lines, texts = [], []

    def draw(n: Node) -> None:
        x0 = left + depth[id(n)] * sx
        if n.children:
            ys = [y[id(c)] for c in n.children]
            lines.append(f"M{x0:.1f},{min(ys):.1f}V{max(ys):.1f}")
            for c in n.children:
                lines.append(f"M{x0:.1f},{y[id(c)]:.1f}H{left + depth[id(c)] * sx:.1f}")
                draw(c)
            if n.name and n.parent is not None:
                texts.append(f'<text x="{x0 - 3:.1f}" y="{y[id(n)] - 3:.1f}" text-anchor="end" '
                             f'class="t-tiny t-muted">{esc(n.name)}</text>')
            return
        x, yy = x0 + 5, y[id(n)]
        slot = slots.get(n.name)
        if slot is not None:
            texts.append(f'<rect x="{x:.1f}" y="{yy - 4.5:.1f}" width="9" height="9" rx="1.5" '
                         f'fill="{slot_colour(slot)}"><title>group {slot + 1}</title></rect>')
        if slots:
            x += 13
        if n.name == "Reference":
            note = f' <tspan class="t-muted" font-weight="400">{esc(ref_name)}</tspan>' if ref_name else ""
            texts.append(f'<text x="{x:.1f}" y="{yy + size * 0.35:.1f}" class="t-ink" font-size="{size}" '
                         f'font-weight="600">Reference{note}</text>')
        else:
            texts.append(f'<text x="{x:.1f}" y="{yy + size * 0.35:.1f}" class="t-ink" font-size="{size}">'
                         f'{esc(n.name)}</text>')

    draw(root)
    parts.append(f'<path d="{" ".join(lines)}" class="branch"/>')
    parts.extend(texts)
    if max_depth > 0:
        bar = _nice_ceil(max_depth / 5)
        y_bar = height - 24
        parts.append(f'<line x1="{left}" y1="{y_bar}" x2="{left + bar * sx:.1f}" y2="{y_bar}" class="branch"/>')
        text = f"{bar:g} substitutions per site"
        if snp_sites:
            count = bar * snp_sites
            shown = f"{round(count):,}" if count >= 10 else f"{count:.2g}"
            text += f" (about {shown} SNP{'s' if count != 1 else ''})"
        parts.append(f'<text x="{left}" y="{y_bar + 13}" class="t-tiny t-muted">{text}</text>')
    parts.append("</svg>")
    return "".join(parts)


# ---------------------------------------------------------------------------------------------------------------
# Figure: the distance heatmap
# ---------------------------------------------------------------------------------------------------------------

def distance_bins(top: int) -> list[tuple[int, int]]:
    """Colour classes for distances 0..top: 0 on its own, then up to four classes whose upper limits are spread
    evenly on a log scale (the quarter powers of the largest distance, rounded down to 1, 2, 3 or 5 x 10^k)."""
    if top <= 0:
        return [(0, 0)]
    edges = sorted({_nice_floor(top ** (k / 4)) for k in (1, 2, 3)} | {top})
    bins, low = [(0, 0)], 1
    for edge in edges:
        if edge >= low:
            bins.append((low, edge))
            low = edge + 1
    return bins


def _bin_classes(count: int) -> list[str]:
    """CSS classes for the bins: q0 for zero, then the blue steps spread over q1..q4."""
    if count <= 1:
        return ["q0"]
    k = count - 1
    return ["q0"] + [f"q{4 if k == 1 else round(1 + 3 * i / (k - 1))}" for i in range(k)]


def heatmap(names: list[str], matrix: dict[str, dict[str, int]], groups: list[list[str]] | None = None) -> str:
    """Distance matrix as an SVG heatmap: binned colours, groups of identical genomes as coloured bands on both
    axes and labelled blocks on the right, the exact distance on hover."""
    slots = group_slots(groups or [])
    n = len(names)
    top_value = max((matrix[a][b] for a in names for b in names), default=0)
    bins = distance_bins(top_value)
    classes = _bin_classes(len(bins))

    def cls(d: int) -> str:
        for (lo, hi), c in zip(bins, classes):
            if lo <= d <= hi:
                return c
        return classes[-1]

    cell = max(10, min(24, 960 // max(n, 1)))
    labelled = n <= MAX_LABELLED_CELLS
    font = min(10.5, cell * 0.55)
    name_px = _text_px(max(names, key=len, default=""), 11)
    band = 14 if slots else 0
    left = 24 + name_px + band
    top = 24 + name_px * 0.87 + band
    right = 130 if slots else 20
    labels = [f"{lo}" if lo == hi else f"{lo}–{hi}" for lo, hi in bins]
    legend_px = 44 + sum(16 + _text_px(label, 11) + 18 for label in labels)
    width = max(left + n * cell + right, left + legend_px + 12)
    height = top + n * cell + 44
    parts = [_svg(width, height, "Pairwise SNP distances")]
    for i, a in enumerate(names):
        x, y = left + i * cell, top + i * cell
        if a in slots:
            colour, title = slot_colour(slots[a]), f"{esc(a)}: group {slots[a] + 1}"
            parts.append(f'<rect x="{left - 14:.0f}" y="{y:.0f}" width="10" height="{cell - 1}" fill="{colour}">'
                         f'<title>{title}</title></rect>')
            parts.append(f'<rect x="{x:.0f}" y="{top - 14:.0f}" width="{cell - 1}" height="10" fill="{colour}">'
                         f'<title>{title}</title></rect>')
        parts.append(f'<text x="{left - band - 6:.0f}" y="{y + cell / 2 + 4:.0f}" text-anchor="end" '
                     f'class="t-small t-ink2">{esc(a)}</text>')
        parts.append(f'<text transform="translate({x + cell / 2 + 3:.0f},{top - band - 6:.0f}) rotate(-60)" '
                     f'class="t-small t-ink2">{esc(a)}</text>')
    by_class: dict[str, list[str]] = {c: [] for c in classes}
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            d = matrix[a][b]
            x, y = left + j * cell, top + i * cell
            title = f"{esc(a)} – {esc(b)}: {d} SNP{'s' if d != 1 else ''}"
            items = by_class[cls(d)]
            items.append(f'<rect x="{x:.0f}" y="{y:.0f}" width="{cell - 1}" height="{cell - 1}"><title>{title}'
                         '</title></rect>')
            if labelled:
                items.append(f'<text x="{x + (cell - 1) / 2:.1f}" y="{y + cell / 2 + font * 0.35:.1f}">{d}</text>')
    for c, items in by_class.items():
        parts.append(f'<g class="{c}" text-anchor="middle" font-size="{font:.1f}">{"".join(items)}</g>')
    i = 0  # Blocks of identical genomes on the right
    while i < n:
        slot = slots.get(names[i])
        j = i
        while j + 1 < n and slots.get(names[j + 1]) == slot:
            j += 1
        if slot is not None:
            size = len(groups[slot]) if groups else j - i + 1
            count = f"{j - i + 1}" if j - i + 1 == size else f"{j - i + 1} of {size}"
            parts.append(f'<rect x="{left + n * cell + 8:.0f}" y="{top + i * cell:.0f}" width="4" '
                         f'height="{(j - i + 1) * cell - 1}" fill="{slot_colour(slot)}"/>')
            parts.append(f'<text x="{left + n * cell + 16:.0f}" y="{top + (i + j + 1) / 2 * cell + 4:.0f}" '
                         f'class="t-small t-ink">group {slot + 1} ({count})</text>')
        i = j + 1
    lx, ly = left, top + n * cell + 14
    parts.append(f'<text x="{lx:.0f}" y="{ly + 10:.0f}" class="t-small t-ink2">SNPs:</text>')
    lx += 44
    for label, c in zip(labels, classes):
        parts.append(f'<g class="{c}"><rect x="{lx:.0f}" y="{ly:.0f}" width="12" height="12"/></g>')
        parts.append(f'<text x="{lx + 16:.0f}" y="{ly + 10:.0f}" class="t-small t-ink2">{label}</text>')
        lx += 16 + _text_px(label, 11) + 18
    parts.append("</svg>")
    return "".join(parts)


# ---------------------------------------------------------------------------------------------------------------
# Figure: the genome map (SNP positions and N bases along the reference)
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class Snp:
    chrom: str
    pos: int
    ref: str
    alt: str
    alt_count: int  # Genomes with an alternate allele
    missing: int  # Genomes without a call


def read_vcf(path: Path) -> tuple[list[Snp], int]:
    """SNP records of a VCF (CHROM, POS, REF, ALT and the genotype counts) and the number of genome columns."""
    snps: list[Snp] = []
    genomes = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("##") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if line.startswith("#"):
                genomes = max(0, len(fields) - 9)
                continue
            if len(fields) < 8 or not fields[1].isdigit():
                continue
            calls = [gt.split(":")[0] for gt in fields[9:]]
            alleles = [re.split(r"[/|]", gt) for gt in calls]
            missing = sum(1 for a in alleles if all(x in (".", "") for x in a))
            alt = sum(1 for a in alleles if any(x.isdigit() and int(x) > 0 for x in a))
            snps.append(Snp(fields[0], int(fields[1]), fields[3], fields[4], alt, missing))
    return snps, genomes


def n_per_bin(assemblies: list[Path], sequences: list[tuple[str, int]], bin_bp: int) -> dict[str, list[int]]:
    """N bases per bin of `bin_bp` along each reference sequence, summed over the assemblies (templated assemblies:
    a record named <sample>_<sequence> follows the coordinates of that reference sequence)."""
    lengths = dict(sequences)
    counts = {name: [0] * (length // bin_bp + 1) for name, length in sequences}
    for path in assemblies:
        prefix = split_extension(path.name)[0] + "_"
        for rec in read_records(path):
            name = rec.name[len(prefix):] if rec.name.startswith(prefix) else rec.name
            if name not in lengths:
                continue
            track = counts[name]
            for i in range(len(track)):
                track[i] += rec.seq.count("N", i * bin_bp, (i + 1) * bin_bp)
    return counts


def genome_map(sequences: list[tuple[str, int]], snps: list[Snp], genomes: int,
               n_tracks: dict[str, list[int]] | None = None, n_bin: int = 1000) -> str:
    """One axis per reference sequence (widths proportional to their lengths), a track of SNP positions and,
    when given, a track of N bases per bin on a log scale."""
    width = FIG_WIDTH
    left, right = 118, width - 20
    longest = max((length for _, length in sequences), default=1) or 1
    scale = (right - left) / longest
    pixel_bp = max(1, math.ceil(longest / (right - left)))
    tick_bp = _nice_ceil(longest / 12)
    block = 30 + 38 + (76 if n_tracks is not None else 0)
    height = 8 + block * len(sequences) + 36
    parts = [_svg(width, height, "SNP positions along the reference")]
    binned: dict[str, dict[int, list[Snp]]] = {}
    for s in snps:
        binned.setdefault(s.chrom, {}).setdefault(s.pos // pixel_bp, []).append(s)
    n_max = max((v for track in (n_tracks or {}).values() for v in track), default=0)
    y = 8
    for name, length in sequences:
        end = left + length * scale
        parts.append(f'<text x="{left}" y="{y + 12}" class="t-small t-ink" font-weight="600">{esc(name)} '
                     f'<tspan class="t-muted" font-weight="400">{length:,} bp</tspan></text>')
        ya = y + 22
        parts.append(f'<line x1="{left}" y1="{ya}" x2="{end:.1f}" y2="{ya}" class="axis"/>')
        tick = 0.0
        while tick <= length:
            x = left + tick * scale
            label = f"{tick / 1e6:g} Mb" if tick_bp >= 1e6 else f"{tick / 1e3:g} kb"
            parts.append(f'<line x1="{x:.1f}" y1="{ya}" x2="{x:.1f}" y2="{ya + 4}" class="axis"/>')
            parts.append(f'<text x="{x:.1f}" y="{ya + 15}" text-anchor="middle" class="t-tiny t-muted">{label}'
                         '</text>')
            tick += tick_bp
        ys = ya + 24
        parts.append(f'<text x="{left - 6}" y="{ys + 12}" text-anchor="end" class="t-small t-muted">SNPs</text>')
        parts.append(f'<line x1="{left}" y1="{ys + 18}" x2="{end:.1f}" y2="{ys + 18}" class="grid"/>')
        for _px, group in sorted(binned.get(name, {}).items()):
            x = left + group[0].pos * scale
            if len(group) == 1:
                s = group[0]
                title = (f"{esc(s.chrom)}:{s.pos:,} {esc(s.ref)}&gt;{esc(s.alt)}: alternate allele in "
                         f"{s.alt_count} of {genomes} genomes; "
                         + (f"{s.missing} missing call{'s' if s.missing != 1 else ''}" if s.missing
                            else "called in every genome"))
                missing = s.missing > 0
            else:
                incomplete = sum(1 for s in group if s.missing)
                title = (f"{len(group)} SNPs at {esc(name)}:{group[0].pos:,}–{group[-1].pos:,}; {incomplete} with "
                         "a missing call. " + " ".join(f"{s.pos:,} {esc(s.ref)}&gt;{esc(s.alt)} ({s.alt_count} alt"
                                                       + (f", {s.missing} missing)" if s.missing else ")")
                                                       for s in group[:8])
                         + (" …" if len(group) > 8 else ""))
                missing = incomplete > 0
            parts.append(f'<line x1="{x:.1f}" y1="{ys}" x2="{x:.1f}" y2="{ys + 18}" '
                         f'class="{"snp-miss" if missing else "snp-all"}"><title>{title}</title></line>')
        y_n = ys + 30
        if n_tracks is not None:
            h = 56
            track = n_tracks.get(name, [])
            unit = f"{n_bin / 1000:g} kb"
            parts.append(f'<text x="{left - 6}" y="{y_n + h / 2:.1f}" text-anchor="end" class="t-small t-muted">'
                         'N bases</text>')
            parts.append(f'<text x="{left - 6}" y="{y_n + h / 2 + 12:.1f}" text-anchor="end" '
                         f'class="t-tiny t-muted">per {unit}, all assemblies</text>')
            parts.append(f'<line x1="{left}" y1="{y_n + h}" x2="{end:.1f}" y2="{y_n + h}" class="axis"/>')
            for i, v in enumerate(track):
                if not v:
                    continue
                hh = h * math.log10(1 + v) / math.log10(1 + max(n_max, 1))
                x = left + i * n_bin * scale
                bar_w = max(1.5, n_bin * scale - 0.5)
                parts.append(f'<rect x="{x:.1f}" y="{y_n + h - hh:.1f}" width="{bar_w:.1f}" '
                             f'height="{hh:.1f}" class="nbar"><title>{esc(name)}:{i * n_bin + 1:,}–'
                             f'{min(length, (i + 1) * n_bin):,}: {v:,} N bases summed over the assemblies</title>'
                             '</rect>')
            if name == sequences[0][0]:
                parts.append(f'<text x="{left + 4}" y="{y_n + 9}" class="t-tiny t-muted">max {n_max:,} N per '
                             f'{unit} (log scale)</text>')
        y += block
    ly = height - 14
    legend = [("snp-all", "SNP called in every genome"), ("snp-miss", "SNP with a missing call in some genomes")]
    lx = left
    for c, label in legend:
        parts.append(f'<line x1="{lx + 4}" y1="{ly - 10}" x2="{lx + 4}" y2="{ly + 2}" class="{c}"/>')
        parts.append(f'<text x="{lx + 12}" y="{ly}" class="t-small t-ink2">{label}</text>')
        lx += 12 + _text_px(label, 11) + 24
    parts.append("</svg>")
    return "".join(parts)


def read_sequences(path: Path) -> list[tuple[str, int]]:
    return [(rec.name, len(rec.seq)) for rec in read_records(path)]


# ---------------------------------------------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------------------------------------------

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
    ref_desc = html.escape(Path(ref.get('file') or s.get('reference') or 'the reference').name)
    if ref.get("length"):
        ref_desc += f", {ref['length']:,} bp"
    if s.get("baiting") == "bbduk":
        hdist = s.get("hdist", 2)  # 2 before 0.3.3
        mismatches = {0: "no mismatch", 1: "up to one mismatch"}.get(hdist, f"up to {hdist} mismatches")
        parts.append(f"Reads sharing a {s.get('kmer', 31)}-mer ({mismatches}) with the reference "
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


# ---------------------------------------------------------------------------------------------------------------
# Page
# ---------------------------------------------------------------------------------------------------------------

TOKENS_LIGHT = """color-scheme:light;--bg:#f9f9f7;--surface:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;
--grid:#e1e0d9;--axis:#c3c2b7;--border:rgba(11,11,11,.10);
--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a;--s4:#eda100;--s5:#e87ba4;--s6:#008300;--s7:#4a3aa7;--s8:#e34948;
--good:#0ca30c;--warning:#fab219;--crit:#d03b3b;
--q0:#e1e0d9;--q1:#86b6ef;--q2:#3987e5;--q3:#1c5cab;--q4:#0d366b;
--qt0:#52514e;--qt1:#0b0b0b;--qt2:#fff;--qt3:#fff;--qt4:#fff;
--bad-bg:#fbe3e3;--bad-fg:#8f1d1d;--warn-bg:#fdf1cf;--warn-fg:#6b4a00;--info-bg:#e1ecfa;--info-fg:#1c4f8f"""
TOKENS_DARK = """color-scheme:dark;--bg:#0d0d0d;--surface:#1a1a19;--ink:#fff;--ink2:#c3c2b7;--muted:#898781;
--grid:#2c2c2a;--axis:#383835;--border:rgba(255,255,255,.10);
--s1:#3987e5;--s2:#d95926;--s3:#199e70;--s4:#c98500;--s5:#d55181;--s6:#008300;--s7:#9085e9;--s8:#e66767;
--good:#0ca30c;--warning:#fab219;--crit:#d03b3b;
--q0:#2c2c2a;--q1:#184f95;--q2:#2a78d6;--q3:#6da7ec;--q4:#b7d3f6;
--qt0:#c3c2b7;--qt1:#fff;--qt2:#fff;--qt3:#0b0b0b;--qt4:#0b0b0b;
--bad-bg:#4a1f1f;--bad-fg:#ffb4ab;--warn-bg:#45371a;--warn-fg:#f5cf82;--info-bg:#1c3050;--info-fg:#a9c8f5"""

CSS = f"""
:root{{{TOKENS_LIGHT}}}
@media (prefers-color-scheme:dark){{:root:not([data-theme="light"]){{{TOKENS_DARK}}}}}
:root[data-theme="dark"]{{{TOKENS_DARK}}}
*{{box-sizing:border-box}}
:root{{print-color-adjust:exact;-webkit-print-color-adjust:exact}}
body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}}
main{{max-width:1180px;margin:0 auto;padding:24px 16px 80px}}
h1{{font-size:26px;line-height:1.2;margin:0 0 6px}}
h1 b{{color:var(--s1)}}
h2{{font-size:21px;margin:44px 0 10px;padding-top:12px;border-top:1px solid var(--grid)}}
h3{{font-size:17px;margin:26px 0 8px}}
p,li{{max-width:900px}}
a{{color:var(--s1)}}
code{{font-size:13px;background:var(--surface);border:1px solid var(--border);border-radius:3px;padding:0 3px}}
.meta{{color:var(--ink2);font-size:14px;margin:0 0 16px}}
.tiles{{display:flex;flex-wrap:wrap;gap:10px;margin:12px 0}}
.tile{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:10px 14px;
min-width:150px;flex:1 1 150px}}
.tile-v{{font-size:26px;font-weight:600;line-height:1.1}}
.tile-l{{font-size:13px;color:var(--ink2);margin-top:4px}}
.tile-s{{font-size:12px;color:var(--muted);margin-top:2px}}
.card{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:14px 16px;margin:14px 0}}
.fig{{display:block;height:auto;background:var(--surface);border:1px solid var(--border);border-radius:8px;
margin:8px 0;break-inside:avoid}}
.figcap{{font-size:13px;color:var(--ink2);max-width:900px;margin:4px 0 16px}}
.figcap b{{color:var(--ink)}}
.two{{display:grid;grid-template-columns:1fr 1fr;gap:16px;align-items:start}}
@media(max-width:900px){{.two{{grid-template-columns:1fr}}}}
.t-small{{font-size:11px;fill:var(--ink2)}}.t-tiny{{font-size:9.5px;fill:var(--ink2)}}
.t-ink{{fill:var(--ink)}}.t-ink2{{fill:var(--ink2)}}.t-muted{{fill:var(--muted)}}
.fig text{{font-family:system-ui,-apple-system,"Segoe UI",sans-serif}}
.grid{{stroke:var(--grid);stroke-width:1}}.axis{{stroke:var(--axis);stroke-width:1}}
.thresh{{stroke:var(--ink2);stroke-width:1}}
.branch{{stroke:var(--ink2);stroke-width:1.5;fill:none;stroke-linecap:square}}
.bar{{fill:var(--s1)}}.bar0{{fill:var(--axis)}}.nbar{{fill:var(--muted)}}
.snp-all{{stroke:var(--s1);stroke-width:1.5}}.snp-miss{{stroke:var(--s2);stroke-width:1.5}}
g.q0 rect{{fill:var(--q0)}}g.q1 rect{{fill:var(--q1)}}g.q2 rect{{fill:var(--q2)}}g.q3 rect{{fill:var(--q3)}}
g.q4 rect{{fill:var(--q4)}}
g.q0 text{{fill:var(--qt0)}}g.q1 text{{fill:var(--qt1)}}g.q2 text{{fill:var(--qt2)}}g.q3 text{{fill:var(--qt3)}}
g.q4 text{{fill:var(--qt4)}}
.fig text{{font-variant-numeric:tabular-nums}}
.tablewrap{{overflow-x:auto;margin:8px 0 16px;border:1px solid var(--border);border-radius:8px;
background:var(--surface)}}
table{{border-collapse:collapse;font-size:13px;width:100%}}
th,td{{padding:5px 8px;text-align:left;vertical-align:top;border-bottom:1px solid var(--grid);white-space:nowrap}}
th{{position:sticky;top:0;background:var(--surface);color:var(--ink2);font-weight:600}}
td{{font-variant-numeric:tabular-nums}}
tbody tr:hover{{background:var(--bg)}}
table.samples{{font-size:12px}}table.samples th,table.samples td{{padding:4px 6px}}
table.samples th{{cursor:pointer;user-select:none}}
table.samples th:after{{content:" \\2195";color:var(--muted);font-size:11px}}
table.samples th:hover{{color:var(--ink)}}
table.samples th.num,table.samples td.num{{text-align:right}}
table.samples td.note{{white-space:normal;min-width:200px}}
.bad{{background:var(--bad-bg);color:var(--bad-fg)}}
.warn{{background:var(--warn-bg);color:var(--warn-fg)}}
.info{{background:var(--info-bg);color:var(--info-fg)}}
.key span{{display:inline-block;padding:1px 8px;border-radius:4px;margin-right:8px;font-size:13px}}
.swatch{{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:6px;vertical-align:middle}}
.tree{{background:#fff;border-radius:8px;padding:8px;overflow-x:auto;border:1px solid var(--border)}}
.tree svg{{max-width:100%;height:auto}}
ul.groups li{{margin-bottom:4px}}
.methods{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:12px 16px}}
dl{{display:grid;grid-template-columns:max-content 1fr;gap:4px 16px;font-size:13px}}
dt{{color:var(--muted)}}
dd{{margin:0;word-break:break-all}}
@media print{{
  body{{background:#fff}} main{{max-width:none;padding:0}}
  table.samples th{{position:static}} .tablewrap{{overflow:visible}}
  h2{{break-after:avoid}} .figcap{{break-before:avoid}}
}}
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


def tile(value: str, label: str, sub: str = "") -> str:
    sub_html = f'<div class="tile-s">{esc(sub)}</div>' if sub else ""
    return (f'<div class="tile"><div class="tile-v">{esc(value)}</div><div class="tile-l">{esc(label)}</div>'
            f"{sub_html}</div>")


class _Figures:
    """Numbered figure captions, in page order."""

    def __init__(self) -> None:
        self.count = 0

    def caption(self, text: str) -> str:
        self.count += 1
        return f'<p class="figcap"><b>Figure {self.count}.</b> {text}</p>'


def _samples_table(rows: list[dict[str, str]]) -> str:
    out = ['<div class="tablewrap"><table class="samples"><thead><tr>'
           + "".join(f'<th class="{"num" if key in NUMERIC else ""}">{esc(label)}</th>'
                     for key, label in TABLE_COLUMNS) + "</tr></thead><tbody>"]
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
    return "".join(out)


def build_report(output: Path) -> str:
    info = json.loads((output / "run_info.json").read_text())
    rows = _read_tsv(output / "summary.tsv") if (output / "summary.tsv").exists() else []
    settings = info.get("settings") or {}
    comparison = info.get("comparison") or {}
    ref = info.get("reference") or {}
    ref_name = Path(str(ref.get("file") or settings.get("reference") or "")).name
    ok = sum(1 for r in rows if r.get("Status") == "ok")
    failed = len(rows) - ok
    flagged = sum(1 for r in rows if r.get("Status") == "ok" and r.get("Note"))
    figures = _Figures()
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

    # The comparison's files, read before the page so that the overview can count the distinct genomes
    distances = located("distances")
    names: list[str] = []
    matrix: dict[str, dict[str, int]] = {}
    order: list[str] = []
    tree: Node | None = None
    in_tree_order = False
    if distances is not None and distances.exists():
        names, matrix = _read_matrix(distances)
        order = names
        tree_file = located("tree")
        if tree_file is not None and tree_file.exists():
            try:
                tree = parse(tree_file.read_text())
                ladderize(tree)
            except Exception:  # noqa: BLE001 - any unreadable tree: the drawn tree.svg is the fallback
                tree = None
            if tree is not None:
                leaves = [leaf.name for leaf in tree.leaves()]
                if sorted(leaves) == sorted(names):
                    order, in_tree_order = leaves, True
    core_snps = comparison.get("core_snps")
    groups = identical_groups(order, matrix) if order and core_snps != 0 else []
    slots = group_slots(groups)
    method = comparison.get("method", "")
    method_name = {"ska": "SKA2", "parsnp": "Parsnp"}.get(method, method)
    tree_tool = {"fasttree": "FastTree", "iqtree": "IQ-TREE"}.get(comparison.get("tree_method", ""), "")

    out = ["<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\">",
           '<meta name="viewport" content="width=device-width,initial-scale=1">',
           '<meta name="color-scheme" content="light dark">',  # Its own dark theme: no forced darkening
           f"<title>BACoN report</title><style>{CSS}</style></head><body><main>",
           "<h1><b>BACoN</b> report</h1>",
           f'<p class="meta">{esc(Path(str(output)).name)} · {esc(str(info.get("started", "")))} · '
           f'BACoN {esc(str(info.get("bacon_version", "")))} · last run {info.get("duration_s", "?")} s</p>']

    # Overview
    tiles = [tile(f"{ok}/{len(rows)}", "samples assembled", f"{failed} failed" if failed else ""),
             tile(str(flagged), "assembled with a note",
                  "low depth, length or N bases" if flagged else ""),
             tile(f"{ref.get('length', 0):,} bp" if ref.get("length") else "?", "reference",
                  f"{ref_name}" + (f", {ref['sequences']} sequence(s)" if ref.get("sequences") else ""))]
    if order:
        tiles.append(tile(str(core_snps if core_snps is not None else "–"), "SNP sites",
                          f"{len(order)} genomes, {distinct_count(order, matrix)} distinct" if core_snps else
                          f"{len(order)} genomes"))
    else:
        tiles.append(tile(str(core_snps if core_snps is not None else "–"), "SNP sites"))
    out.append(f'<div class="tiles">{"".join(tiles)}</div>')
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
    out.append(_samples_table(rows))
    if rows:
        depth_chart = bar_chart(depth_bars(rows), "Depth after filtering, per sample", unit="x",
                                threshold=LOW_DEPTH, threshold_text=f"{LOW_DEPTH}x flag")
        n_chart = bar_chart(n_bars(rows), "N bases per assembly", integers=True)
        out.append(f'<div class="two"><div>{depth_chart}'
                   + figures.caption("Estimated depth of the filtered reads over the reference, per sample. "
                                     f"The line is the {LOW_DEPTH}x threshold of the table's depth flag; samples "
                                     "below it are labelled. Failed samples are listed without a bar. Hover a bar "
                                     "for the read counts and the note.")
                   + f"</div><div>{n_chart}"
                   + figures.caption("N bases in each assembly, largest first (the table flags every assembly "
                                     f"with N bases; the {N_LABELS} largest counts are labelled). Hover a bar for "
                                     "the fraction of the assembly.")
                   + "</div></div>")

    # Comparison
    if order:
        data = "Parsnp core-genome alignment" if method == "parsnp" else f"{method_name} SNPs"
        svg = located("tree_svg")
        out.append("<h2>Tree</h2>")
        if tree is not None:
            sites = core_snps if method != "parsnp" and isinstance(core_snps, int) else None
            out.append(tree_svg(tree, slots, ref_name, sites))
            squares = (" Squares mark the genomes with no SNP between them (one colour per group, as in the "
                       "heatmap below)." if groups else "")
            bar_text = (" The scale bar is in substitutions per SNP site, with the equivalent number of SNPs."
                        if sites else " The scale bar is in substitutions per site.")
            out.append(figures.caption(f"{esc(data)}; {esc(tree_tool)}, midpoint-rooted and ladderized. Numbers "
                                       f"on the internal branches are supports.{squares}{bar_text}"))
        elif svg is not None and svg.exists():
            out.append(f'<div class="tree">{svg.read_text()}</div>')
            out.append(figures.caption(f"{esc(data)}; {esc(tree_tool)}, midpoint-rooted; internal labels are "
                                       "supports."))
        elif core_snps == 0:
            why = ("no SNP site is shared by all the genomes; see --ska-min-freq"
                   if method == "ska" else "the alignment has no SNP site")
            out.append(f"<p>No tree: {why}.</p>")
        else:
            out.append("<p>No tree.</p>")
        out.append("<h2>SNP distances</h2>")
        order_text = "in tree order" if in_tree_order else "in the order of the distance table"
        if len(order) <= MAX_HEATMAP_GENOMES:
            out.append(heatmap(order, matrix, groups))
            bands = (" Coloured bands on both axes and the blocks on the right mark the groups of identical "
                     "genomes." if groups else "")
            out.append(figures.caption(f"Pairwise SNP distances, {order_text}. Colour classes are spread on a "
                                       "log scale over the range of the distances (legend); hover a cell for the "
                                       f"exact distance of its pair.{bands}"))
        else:
            shown = distances.relative_to(output) if distances.is_relative_to(output) else distances
            out.append(f"<p>{len(order)} genomes: too many for a heatmap in this page; the distances are in "
                       f"<code>{esc(str(shown))}</code>.</p>")
        if core_snps == 0:
            out.append("<p>No SNP site was compared: the distances say nothing about identity.</p>")
        else:
            if groups:
                out.append('<h3>Identical genomes</h3><p class="meta">No SNP between any two genomes of a group '
                           '(positions with N or a gap are not compared).</p><ul class="groups">' + "".join(
                               f'<li><span class="swatch" style="background:{slot_colour(i)}"></span>'
                               f"<b>group {i + 1}</b> ({len(g)}): {esc(', '.join(g))}</li>"
                               for i, g in enumerate(groups)) + "</ul>")
                if len(groups) > GROUP_COLOURS:
                    out.append(f'<p class="meta">Groups beyond the {GROUP_COLOURS}th share the grey colour.</p>')
            out.append(f"<p>{len(order)} genomes, {distinct_count(order, matrix)} distinct at the SNP sites "
                       "compared.</p>")

        # Genome map
        vcf = located("vcf")
        reference = output / "reference.fasta"
        if vcf is not None and vcf.exists() and reference.exists():
            out.append(_genome_map_section(output, vcf, reference, settings, figures, len(order)))
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


def _genome_map_section(output: Path, vcf: Path, reference: Path, settings: dict, figures: _Figures,
                        genomes: int) -> str:
    try:
        sequences = read_sequences(reference)
        snps, columns = read_vcf(vcf)
    except Exception as exc:  # noqa: BLE001 - a malformed file loses its figure, not the report
        return f'<p class="meta">No genome map: {esc(str(exc))}.</p>'
    if not sequences:
        return ""
    drawn = sequences
    skipped = ""
    if len(sequences) > MAX_MAP_SEQUENCES:
        longest = set(sorted(sequences, key=lambda s: -s[1])[:MAX_MAP_SEQUENCES])
        drawn = [s for s in sequences if s in longest]
        names = {name for name, _ in drawn}
        left_out = sum(1 for s in snps if s.chrom not in names)
        skipped = (f" The {len(sequences) - len(drawn)} shorter sequences of the reference ({left_out} SNPs) are "
                   "not drawn.")
    assembler = settings.get("assembler")
    n_tracks = None
    n_bin = 1000
    n_text = ""
    if assembler == "samtools":
        assemblies = sorted((output / "3_assembled" / "all_assemblies").glob("*.fasta"))
        if assemblies:
            total = max(length for _, length in drawn)
            n_bin = int(_nice_ceil(total / 1500)) if total > 1_500_000 else 1000
            try:
                n_tracks = n_per_bin(assemblies, drawn, n_bin)
            except Exception:  # noqa: BLE001 - an unreadable assembly loses the N track only
                n_tracks = None
        if n_tracks is None:
            n_text = " No N track: the assemblies were not found."
        elif not any(any(track) for track in n_tracks.values()):
            n_tracks = None
            n_text = f" No N track: none of the {len(assemblies)} templated assemblies has an N base."
        else:
            n_text = (f" The N track sums the N bases of the {len(assemblies)} templated assemblies per "
                      f"{n_bin / 1000:g} kb (log scale), where the consensus follows the reference coordinates.")
    else:
        n_text = (f" No N track: {assembler or 'the'} assemblies are de novo, so their coordinates do not follow "
                  "the reference.")
    missing = sum(1 for s in snps if s.missing)
    pixel_bp = max(1, math.ceil(max(length for _, length in drawn) / (FIG_WIDTH - 138)))
    binned = (f" Ticks closer than {pixel_bp:,} bp are merged; hover shows the SNPs of a tick."
              if pixel_bp > 1 else "")
    genomes_text = f"{columns} genomes" if columns else f"{genomes} genomes"
    return ("<h2>Genome map</h2>" + genome_map(drawn, snps, columns or genomes, n_tracks, n_bin)
            + figures.caption(f"SNP positions along the reference ({len(snps):,} records of the VCF, {missing:,} "
                              f"with a missing call in at least one of the {genomes_text}). Hover a tick for the "
                              "position, the alleles and the number of genomes with the alternate allele."
                              f"{binned}{n_text}{skipped}"))


def write_report(output: Path) -> Path:
    path = output / "report.html"
    tmp = path.with_suffix(".html.tmp")
    tmp.write_text(build_report(output), encoding="utf-8")
    tmp.replace(path)
    return path


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python -m bacon.report OUTPUT_FOLDER")
    if not (Path(sys.argv[1]) / "run_info.json").is_file():
        sys.exit(f"{sys.argv[1]}: no run_info.json (not a BACoN output folder)")
    print(write_report(Path(sys.argv[1])))
