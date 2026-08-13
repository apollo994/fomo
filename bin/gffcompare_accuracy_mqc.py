#!/usr/bin/env python3
"""
Build a MultiQC custom-content scatter plot of gffcompare Sensitivity vs
Precision, replacing MultiQC's native gffcompare "Accuracy values" plot.

The native plot hard-codes each scatter point's name to the feature *level*
(Base/Exon/...), so hovering a dot never reveals which sample it belongs to.
This builds the same plot — six switchable level datasets, Precision vs
Sensitivity on 0..1 axes — but encodes:

  * marker SYMBOL  = feature class  (lncRNA → circle, mRNA → square,
                                     decoy → triangle-up)
  * marker COLOUR  = source         (one colour per source species, plus one per
                                     aggregate model: allModels / top3, raw and
                                     collapsed)
  * point NAME     = full sample id (so the tooltip is informative)

Input: one or more gffcompare `.stats` files named
       <target>.from_<source>.<feature_type>[.decoy].gffcompare.stats
Output: a single <prefix>_mqc.json MultiQC custom-content file.
"""
import argparse
import json
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

# Feature levels, in the order gffcompare reports them and the native plot
# lists them in its dropdown.
LEVELS = ["Base", "Exon", "Intron", "Intron_chain", "Transcript", "Locus"]

# Marker symbol per feature class (Plotly symbol names).
SYMBOL_BY_CLASS = {
    "lncRNA": "circle",
    "mRNA":   "square",
    "decoy":  "triangle-up",
}

# Qualitative palette (D3 category10) assigned per source.
PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]

FNAME_RE = re.compile(
    r"^(?P<target>.+?)\.from_(?P<source>.+?)\.(?P<ft>lnc_RNA|mRNA)(?P<decoy>\.decoy)?$"
)

# Aggregate pseudo-sources, in the order they should take palette colours. Listing
# them first keeps the four headline models on stable colours across runs with
# different species. Kept in sync with AGGREGATE_IDS in bin/select_top_sources.py
# and subworkflows/local/consensus_top.nf.
AGGREGATE_IDS = ("allModels_raw", "allModels_collapsed", "top3_raw", "top3_collapsed")


def parse_identity(stats_path: str) -> Optional[Tuple[str, str, bool]]:
    """Derive (source, feature_type, decoy) from a gffcompare stats filename."""
    base = os.path.basename(stats_path)
    base = re.sub(r"\.gffcompare\.stats$", "", base)
    base = re.sub(r"\.stats$", "", base)            # tolerate either suffix
    base = re.sub(r"\.gffcompare$", "", base)
    m = FNAME_RE.match(base)
    if not m:
        return None
    return m.group("source"), m.group("ft"), bool(m.group("decoy"))


def feature_class(ft: str, decoy: bool) -> str:
    if decoy:
        return "decoy"
    return "lncRNA" if ft == "lnc_RNA" else "mRNA"


def parse_accuracy(stats_path: str) -> Dict[str, Tuple[float, float]]:
    """Return {level: (sensitivity, precision)} read from a stats file.
    Mirrors MultiQC's native parser: drop '|', normalise 'Intron chain',
    then split — token[0]=level, token[2]=sensitivity, token[3]=precision."""
    out: Dict[str, Tuple[float, float]] = {}
    with open(stats_path, encoding="utf-8") as fh:
        for line in fh:
            if "level:" not in line:
                continue
            toks = line.replace("|", "").replace("Intron chain", "Intron_chain").split()
            if len(toks) < 4:
                continue
            level = toks[0]
            try:
                sens = float(toks[2])
                prec = float(toks[3])
            except ValueError:
                continue
            out[level] = (sens, prec)
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stats", nargs="+", help="gffcompare .stats files")
    p.add_argument("--output", required=True, metavar="JSON",
                   help="output MultiQC custom-content JSON (name it *_mqc.json)")
    p.add_argument("--section-name", default="Accuracy values")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    # Collect (source, ft, decoy, accuracy) per file
    records = []
    for path in args.stats:
        ident = parse_identity(path)
        if ident is None:
            print(f"WARNING: cannot parse identity from '{path}', skipping", file=sys.stderr)
            continue
        source, ft, decoy = ident
        acc = parse_accuracy(path)
        if not acc:
            print(f"WARNING: no accuracy lines in '{path}', skipping", file=sys.stderr)
            continue
        records.append((source, ft, decoy, acc))

    if not records:
        print("ERROR: no usable gffcompare stats files", file=sys.stderr)
        return 1

    # Stable colour assignment: the aggregate pseudo-sources first, in a fixed
    # order, then the real species alphabetically — so the four headline models keep
    # the same colours no matter which species a run happens to include.
    present = {r[0] for r in records}
    sources = [s for s in AGGREGATE_IDS if s in present] + sorted(present - set(AGGREGATE_IDS))
    colour_by_source = {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(sources)}

    # One data dict per level (aligned with data_labels) → multi-dataset plot
    data: List[Dict[str, dict]] = []
    for level in LEVELS:
        level_points: Dict[str, dict] = {}
        for source, ft, decoy, acc in records:
            if level not in acc:
                continue
            sens, prec = acc[level]
            name = f"{source}.{ft}" + (".decoy" if decoy else "")
            level_points[name] = {
                "x": round(sens / 100.0, 4),
                "y": round(prec / 100.0, 4),
                "color": colour_by_source[source],
                "marker_symbol": SYMBOL_BY_CLASS[feature_class(ft, decoy)],
                "name": name,
            }
        data.append(level_points)

    doc = {
        "id": "gffcompare_accuracy_custom",
        "section_name": args.section_name,
        "description": (
            "Sensitivity vs Precision per feature level (use the dropdown to "
            "switch level). Marker symbol = feature class "
            "(lncRNA ●, mRNA ■, decoy ▲); colour = source species "
            "or aggregate model (allModels / top3, each raw and collapsed). "
            "Hover a point for its sample name."
        ),
        "plot_type": "scatter",
        "pconfig": {
            "id": "gffcompare_accuracy_custom_plot",
            "title": "GffCompare: Accuracy values",
            "xlab": "Sensitivity",
            "ylab": "Precision",
            "xmin": 0, "xmax": 1, "ymin": 0, "ymax": 1,
            "data_labels": [
                {"name": lvl.replace("_", " "), "xlab": "Sensitivity", "ylab": "Precision"}
                for lvl in LEVELS
            ],
        },
        "data": data,
    }

    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2)

    n_pts = sum(len(d) for d in data)
    print(f"gffcompare_accuracy_mqc: {len(records)} samples × {len(LEVELS)} levels "
          f"= {n_pts} points; {len(sources)} colours", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
