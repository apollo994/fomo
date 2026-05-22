#!/usr/bin/env python3
"""
Plot distributions of NM, AS, de from alignment_metrics.tsv.

Layout:
    - 3 columns: NM, AS, de
    - 1 row per (source, target) pair
    - Each subplot: per-genetype (mrna, lnc, decoy) violin distribution

Usage:
    03_plot_alignment_metrics.py alignment_metrics.tsv -o alignment_metrics.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

METRICS = ["NM", "AS", "de"]
GENETYPE_ORDER = ["mrna", "lnc", "decoy"]
GENETYPE_COLOR = {"mrna": "#1f77b4", "lnc": "#2ca02c", "decoy": "#d62728"}


def violin(ax, data_by_group, order, colors):
    positions = list(range(len(order)))
    series = [np.asarray(data_by_group.get(g, []), dtype=float) for g in order]
    series = [s[~np.isnan(s)] if s.size else s for s in series]

    valid_idx = [i for i, s in enumerate(series) if s.size > 0]
    if valid_idx:
        parts = ax.violinplot(
            [series[i] for i in valid_idx],
            positions=[positions[i] for i in valid_idx],
            showmeans=False,
            showmedians=True,
            widths=0.8,
        )
        for i, body in zip(valid_idx, parts["bodies"]):
            body.set_facecolor(colors[order[i]])
            body.set_edgecolor("black")
            body.set_alpha(0.6)
        for key in ("cbars", "cmins", "cmaxes", "cmedians"):
            if key in parts:
                parts[key].set_color("black")
                parts[key].set_linewidth(0.8)

    ax.set_xticks(positions)
    ax.set_xticklabels(order, fontsize=8)
    ymax = max((s.max() for s in series if s.size), default=1.0)
    for i, s in enumerate(series):
        ax.text(positions[i], ymax, f"n={s.size}", ha="center", va="bottom", fontsize=6)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", type=Path, help="alignment_metrics.tsv from extract script")
    ap.add_argument("-o", "--output", type=Path, default=Path("alignment_metrics.png"))
    ap.add_argument("--row-height", type=float, default=2.2)
    ap.add_argument("--col-width", type=float, default=3.2)
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    for c in METRICS:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # group rows by target so all comparisons sharing a target are adjacent
    pairs = sorted(
        df[["source", "target"]].drop_duplicates().itertuples(index=False, name=None),
        key=lambda st: (st[1], st[0]),
    )
    n_rows = len(pairs)
    n_cols = len(METRICS)

    fig = plt.figure(figsize=(args.col_width * n_cols, args.row_height * n_rows))
    subfigs = np.atleast_1d(fig.subfigures(n_rows, 1, hspace=0.4))

    for r, (source, target) in enumerate(pairs):
        sub = df[(df["source"] == source) & (df["target"] == target)]
        sf = subfigs[r]
        sf.suptitle(
            f"target: {target}   ←   source: {source}",
            fontsize=11, fontweight="bold",
        )
        axs = sf.subplots(1, n_cols)
        for c, metric in enumerate(METRICS):
            ax = axs[c]
            data_by_group = {g: sub.loc[sub["genetype"] == g, metric].dropna().values
                             for g in GENETYPE_ORDER}
            violin(ax, data_by_group, GENETYPE_ORDER, GENETYPE_COLOR)
            ax.set_title(metric, fontsize=10, fontweight="bold")
            ax.tick_params(axis="y", labelsize=7)
            ax.grid(axis="y", linestyle=":", alpha=0.4)

    fig.suptitle("Minimap transfer alignment metrics by target ← source and genetype",
                 fontsize=12, fontweight="bold")
    fig.savefig(args.output, bbox_inches="tight", dpi=args.dpi)
    print(f"wrote {args.output} ({n_rows} rows × {n_cols} cols, {len(df)} total rows)")


if __name__ == "__main__":
    main()
