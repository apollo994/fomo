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
import sys
from pathlib import Path
from typing import Dict, List

# Nextflow only puts bin/ on PATH, not on Python's import path — but it stages the
# whole directory, so a __file__-relative insert finds the sibling module. The
# identity grammar, AGGREGATE_IDS, the level list and the stats parser live there
# because this script and select_top_sources.py must agree on all of them.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fomo_stats import (  # noqa: E402
    AGGREGATE_IDS,
    LEVELS,
    feature_class,
    parse_accuracy,
    parse_identity,
)

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

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stats", nargs="+", help="gffcompare .stats files")
    p.add_argument("--output", required=True, metavar="JSON",
                   help="output MultiQC custom-content JSON (name it *_mqc.json)")
    p.add_argument("--section-name", default="Accuracy values")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    # Collect (identity, accuracy) per file
    records = []
    for path in args.stats:
        ident = parse_identity(path)
        if ident is None:
            print(f"WARNING: cannot parse identity from '{path}', skipping", file=sys.stderr)
            continue
        acc = parse_accuracy(path)
        if not acc:
            print(f"WARNING: no accuracy lines in '{path}', skipping", file=sys.stderr)
            continue
        records.append((ident, acc))

    if not records:
        # Every comparison for this target had no accuracy lines — happens when the
        # target's reference GFF3 is non-empty pre-filter but filters down to zero
        # spliced transcripts (FILTER_TARGET keeps multi-exon only), so gffcompare
        # falls back to combine-mode-style output with no Sensitivity/Precision
        # section for any source. That is a real "nothing to plot" outcome for this
        # target, not a pipeline failure — fall through and emit a validly-shaped but
        # empty scatter (every `data` entry ends up `{}`, same as the per-level empty
        # case below when `records` is non-empty but no source has a given level).
        print("WARNING: no usable gffcompare stats files (every comparison had an "
              "empty or feature-less reference after filtering) — emitting an empty "
              "accuracy plot", file=sys.stderr)

    # Stable colour assignment: the aggregate pseudo-sources first, in the fixed
    # AGGREGATE_IDS order, then the real species alphabetically — so the four headline
    # models keep the same colours no matter which species a run happens to include.
    present = {ident.source for ident, _acc in records}
    sources = [s for s in AGGREGATE_IDS if s in present] + sorted(present - set(AGGREGATE_IDS))
    colour_by_source = {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(sources)}

    # One data dict per level (aligned with data_labels) → multi-dataset plot
    data: List[Dict[str, dict]] = []
    for level in LEVELS:
        level_points: Dict[str, dict] = {}
        for ident, acc in records:
            if level not in acc:
                continue
            sens, prec = acc[level]
            name = f"{ident.source}.{ident.feature_type}" + (".decoy" if ident.decoy else "")
            level_points[name] = {
                "x": round(sens / 100.0, 4),
                "y": round(prec / 100.0, 4),
                "color": colour_by_source[ident.source],
                "marker_symbol": SYMBOL_BY_CLASS[
                    feature_class(ident.feature_type, ident.decoy)
                ],
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
