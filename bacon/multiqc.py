"""MultiQC custom-content files (https://docs.seqera.io/multiqc/custom_content).

Three files in the output folder, found by `multiqc OUTPUT_FOLDER`:
  bacon_samples_mqc.json    table: reads, depth, assembly and status of each sample
  bacon_reads_mqc.json      bar graph: each sample's bases kept, baited but filtered out, and off-target
  bacon_distances_mqc.json  heatmap: pairwise SNP distances (when the samples were compared)
"""

from __future__ import annotations

import json
from pathlib import Path

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


def sample_table(rows: list[dict[str, str]]) -> dict:
    data = {}
    for r in rows:
        values = {k: r.get(k, "") if k in ("Status", "Note") else _number(r.get(k, "")) for k in SAMPLE_HEADERS}
        data[r["Sample"]] = {k: v for k, v in values.items() if v not in (None, "")}
    return {
        "id": "bacon_samples",
        "section_name": "BACoN: samples",
        "description": "Reads baited and kept, depth and assembly of each sample (BACoN summary.tsv).",
        "plot_type": "table",
        "pconfig": {"id": "bacon_samples_table", "title": "BACoN: samples", "namespace": "BACoN"},
        "headers": SAMPLE_HEADERS,
        "data": data,
    }


def reads_bargraph(rows: list[dict[str, str]]) -> dict | None:
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
        "id": "bacon_reads",
        "section_name": "BACoN: bases",
        "description": "Bases of each sample: kept for the assembly, matching the reference but filtered out "
                       "(short, low quality, or above the target depth), and not matching the reference.",
        "plot_type": "bargraph",
        "categories": {"Kept": {"color": "#2f7ebc"}, "Baited, filtered out": {"color": "#9ecae1"},
                       "Baited": {"color": "#2f7ebc"}, "Off-target": {"color": "#d9d9d9"}},
        "pconfig": {"id": "bacon_reads_plot", "title": "BACoN: bases", "ylab": "Bases",
                    "cpswitch_counts_label": "Bases"},
        "data": data,
    }


def distance_heatmap(path: Path, tree: Path | None = None) -> dict:
    lines = path.read_text().splitlines()
    names = lines[0].split("\t")[1:]
    rows = {line.split("\t")[0]: [int(x) for x in line.split("\t")[1:]] for line in lines[1:]}
    order = names
    if tree is not None and tree.exists():  # Rows and columns in tree order, like report.html
        leaves = [leaf.name for leaf in parse(tree.read_text()).leaves()]
        if sorted(leaves) == sorted(names):
            order = leaves
    index = {n: i for i, n in enumerate(names)}
    matrix = [[rows[a][index[b]] for b in order] for a in order]
    return {
        "id": "bacon_distances",
        "section_name": "BACoN: SNP distances",
        "description": "Pairwise SNP distances between the assemblies and the reference (clustered view first; "
                       "switch to sorted by sample above the plot).",
        "plot_type": "heatmap",
        "pconfig": {"id": "bacon_distances_heatmap", "title": "BACoN: SNP distances", "square": True, "min": 0,
                    "cluster_switch_clustered_active": True,
                    "colstops": [[0, "#ffffd9"], [0.25, "#a1dab4"], [0.5, "#41b6c4"], [0.75, "#225ea8"],
                                 [1, "#081d58"]]},
        "xcats": order,
        "ycats": order,
        "data": matrix,
    }


def write_multiqc(output: Path, rows: list[dict[str, str]], distances: Path | None,
                  tree: Path | None = None) -> list[Path]:
    """Write the MultiQC files; remove a stale distance file when there is no comparison."""
    written = []
    sections = [("bacon_samples_mqc.json", sample_table(rows)), ("bacon_reads_mqc.json", reads_bargraph(rows)),
                ("bacon_distances_mqc.json", distance_heatmap(distances, tree) if distances and distances.exists()
                 else None)]
    for name, content in sections:
        path = output / name
        if content is None:
            path.unlink(missing_ok=True)
            continue
        path.write_text(json.dumps(content, indent=1) + "\n")
        written.append(path)
    return written
