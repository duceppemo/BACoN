"""MultiQC custom-content files (https://docs.seqera.io/multiqc/custom_content).

Three files in the output folder, found by `multiqc OUTPUT_FOLDER`. Their section ids include the name of the
output folder, so that the files of several runs in one MultiQC search path make separate sections instead of
being merged (MultiQC merges sections with the same id, and keeps only one heatmap).
  bacon_samples_mqc.json    table: reads, depth, assembly and status of each sample, with its metadata
  bacon_reads_mqc.json      bar graph: each sample's bases kept, baited but filtered out, and off-target
  bacon_distances_mqc.json  heatmap: pairwise SNP distances (when the samples were compared)
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

from bacon import BaconError
from bacon.metadata import Metadata, shown_name
from bacon.newick import parse

SAMPLE_HEADERS = {
    "Status": {"title": "Status", "description": "ok, or the step at which the sample failed"},
    "Raw_reads": {"title": "Raw reads", "description": "Input reads", "format": "{:,.0f}", "scale": "Greys",
                  "hidden": True},
    "Baited_reads": {"title": "Baited reads", "description": "Reads matching the reference", "format": "{:,.0f}",
                     "scale": "Blues"},
    "Baited_pct": {"title": "Baited %", "description": "Share of the input bases matching the reference",
                   "suffix": "%", "min": 0, "max": 100, "format": "{:,.1f}", "scale": "Purples"},
    "Filtered_reads": {"title": "Filtered reads", "description": "Reads kept by Filtlong", "format": "{:,.0f}",
                       "scale": "Blues", "hidden": True},
    "Filtered_N50": {"title": "Read N50", "description": "N50 of the filtered reads", "suffix": " bp",
                     "format": "{:,.0f}", "scale": "Greens"},
    "Est_depth": {"title": "Depth", "description": "Filtered bases / genome size", "suffix": "x",
                  "format": "{:,.1f}", "scale": "RdYlGn", "min": 0},
    "Contigs": {"title": "Contigs", "description": "Contigs in the assembly", "format": "{:,.0f}",
                "scale": "Oranges"},
    "Circular_contigs": {"title": "Circular", "description": "Contigs reported as circular (de novo assemblers)",
                         "format": "{:,.0f}"},
    "Assembly_length": {"title": "Length", "description": "Assembly length", "suffix": " bp",
                        "format": "{:,.0f}", "scale": "Greys"},
    "Length_vs_reference": {"title": "x ref", "description": "Assembly length / reference length",
                            "format": "{:,.3f}", "scale": "RdBu"},
    "N_bases": {"title": "N bases", "description": "Templated assembly: bases called N (low or ambiguous support)",
                "format": "{:,.0f}", "scale": "Reds"},
    "Note": {"title": "Note", "description": "Warnings: low depth, unexpected length, N bases, failure reason"},
}


def _number(value: str) -> float | int | str | None:
    if value in ("", "NA"):
        return None
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def _slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", name).strip("_")


def _run_id(run: str) -> str:
    """The output folder's name in the section ids: as it is when it is made of letters, digits and underscores,
    else made so with a short hash of the name added, so that folders whose names differ only in punctuation
    (run-1, run_1, run 1) keep separate sections."""
    slug = _slug(run)
    if slug == run:
        return run or "run"
    return f"{slug or 'run'}_{hashlib.sha256(run.encode()).hexdigest()[:6]}"


def sample_table(rows: list[dict[str, str]], run: str = "", metadata: Metadata | None = None) -> dict:
    """The samples table; the metadata columns (as text) come first, when there are any."""
    headers: dict[str, dict] = {}
    columns = {}
    for column in (metadata.columns if metadata else []):
        key = base = f"meta_{_slug(column).lower() or 'run'}"
        suffix = 1
        while key in headers:  # Two names differing only in punctuation or case (or like a suffixed key)
            suffix += 1
            key = f"{base}_{suffix}"
        title = shown_name(column, ["Sample", *SAMPLE_HEADERS, *(h["title"] for h in SAMPLE_HEADERS.values())])
        headers[key] = {"title": title, "description": f"Metadata: {column}"}
        columns[key] = column
    headers.update(SAMPLE_HEADERS)
    data = {}
    for r in rows:
        values = {k: r.get(k, "") if k in ("Status", "Note") else _number(r.get(k, "")) for k in SAMPLE_HEADERS}
        if metadata is not None:
            values.update({k: metadata.value(r["Sample"], c) for k, c in columns.items()})
        data[r["Sample"]] = {k: v for k, v in values.items() if v not in (None, "")}
    return {
        "id": f"bacon_samples_{_run_id(run)}",
        "section_name": f"BACoN {run}: samples".replace("  ", " "),
        "description": "Reads baited and kept, depth and assembly of each sample (BACoN summary.tsv)"
                       + (", with the sample metadata." if metadata else "."),
        "plot_type": "table",
        "pconfig": {"id": f"bacon_samples_table_{_run_id(run)}", "title": f"BACoN {run}: samples",
                    "namespace": f"BACoN {run}".strip()},
        "headers": headers,
        "data": data,
    }


def reads_bargraph(rows: list[dict[str, str]], run: str = "") -> dict | None:
    data = {}
    for r in rows:
        raw, baited, kept = (_number(r.get(k, "")) for k in ("Raw_bases", "Baited_bases", "Filtered_bases"))
        if not isinstance(baited, int):
            continue
        entry = {}
        if isinstance(kept, int):
            entry["Kept"] = kept
            entry["Baited, filtered out"] = max(0, baited - kept)
        else:
            entry["Baited"] = baited
        if isinstance(raw, int):
            entry["Off-target"] = max(0, raw - baited)
        data[r["Sample"]] = entry
    if not data:
        return None
    return {
        "id": f"bacon_reads_{_run_id(run)}",
        "section_name": f"BACoN {run}: bases".replace("  ", " "),
        "description": "Bases of each sample: kept for the assembly, matching the reference but filtered out "
                       "(short, low quality, or above the target depth), and not matching the reference.",
        "plot_type": "bargraph",
        "categories": {"Kept": {"color": "#2f7ebc"}, "Baited, filtered out": {"color": "#9ecae1"},
                       "Baited": {"color": "#2f7ebc"}, "Off-target": {"color": "#d9d9d9"}},
        "pconfig": {"id": f"bacon_reads_plot_{_run_id(run)}", "title": f"BACoN {run}: bases", "ylab": "Bases",
                    "cpswitch_counts_label": "Bases"},
        "data": data,
    }


def distance_heatmap(path: Path, tree: Path | None = None, run: str = "") -> dict:
    """The distance heatmap, rows and columns in the order of the tree's leaves (like report.html), or of the
    distance table when the tree is missing, unreadable or of other genomes. The distances are shown as integers
    (tt_decimals) by recent MultiQC versions; older ones such as 1.19 ignore it and show two decimals (their
    decimalPlaces is deprecated in recent versions, which warn about it)."""
    lines = path.read_text().splitlines()
    names = lines[0].split("\t")[1:]
    rows = {line.split("\t")[0]: [int(x) for x in line.split("\t")[1:]] for line in lines[1:]}
    order, in_tree_order = names, False
    leaves: list[str] = []
    if tree is not None and tree.exists():
        try:
            leaves = [leaf.name for leaf in parse(tree.read_text()).leaves()]
        except (BaconError, OSError, UnicodeDecodeError, RecursionError):
            leaves = []
        if leaves and sorted(leaves) == sorted(names):
            order, in_tree_order = leaves, True
    index = {n: i for i, n in enumerate(names)}
    matrix = [[rows[a][index[b]] for b in order] for a in order]
    return {
        "id": f"bacon_distances_{_run_id(run)}",
        "section_name": f"BACoN {run}: SNP distances".replace("  ", " "),
        "description": "Pairwise SNP distances between the assemblies and the reference, rows and columns in "
                       + ("tree order" if in_tree_order else "the order of the distance table")
                       + " (recent MultiQC versions also offer a clustered view).",
        "plot_type": "heatmap",
        "pconfig": {"id": f"bacon_distances_heatmap_{_run_id(run)}", "title": f"BACoN {run}: SNP distances",
                    "square": True, "min": 0, "tt_decimals": 0,
                    "colstops": [[0, "#ffffd9"], [0.25, "#a1dab4"], [0.5, "#41b6c4"], [0.75, "#225ea8"],
                                 [1, "#081d58"]]},
        "xcats": order,
        "ycats": order,
        "data": matrix,
    }


def write_multiqc(output: Path, rows: list[dict[str, str]], distances: Path | None,
                  tree: Path | None = None, metadata: Metadata | None = None) -> list[Path]:
    """Write the MultiQC files; remove a stale distance file when there is no comparison."""
    written = []
    run = output.resolve().name
    sections = [("bacon_samples_mqc.json", sample_table(rows, run, metadata)),
                ("bacon_reads_mqc.json", reads_bargraph(rows, run)),
                ("bacon_distances_mqc.json", distance_heatmap(distances, tree, run)
                 if distances and distances.exists() else None)]
    for name, content in sections:
        path = output / name
        if content is None:
            path.unlink(missing_ok=True)
            continue
        tmp = path.with_name(path.name + ".tmp")
        tmp.write_text(json.dumps(content, indent=1) + "\n")
        tmp.replace(path)
        written.append(path)
    return written
