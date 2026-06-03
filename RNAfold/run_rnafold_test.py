#!/usr/bin/env python
"""
Standalone RNAfold sliding-window MFE test.

For a single species (Lycaena_alciphron_282377) it samples 50 sequences from each of
four classes -- lncRNA real, lncRNA decoy, mRNA real, mRNA decoy -- slides windows of
30/60/90 nt (step 10) over each sequence, computes per-window minimum free energy (MFE)
via the ViennaRNA Python API, and writes tidy CSVs + comparison plots.

This is an exploratory experiment, NOT part of the FOMO Nextflow pipeline.

Run from the repo root:
    conda activate rnafold_test
    python RNAfold/run_rnafold_test.py
"""

import random
from pathlib import Path

import pandas as pd
import matplotlib

matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

import helpers

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
SEED = 42
N_PER_CLASS = 50
WINDOW_SIZES = [30, 60, 90, 120]
STEP = 3
# STEP = 10

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
FASTA_DIR = REPO_ROOT / "results" / "rename_fasta_headers"

SUBSET_DIR = HERE / "subsets"
RESULTS_DIR = HERE / "results"
PLOTS_DIR = HERE / "plots"

SPECIES = "Lycaena_alciphron_282377"

# (rna_class, source, fasta filename)
INPUTS = [
    ("lncRNA", "real", f"{SPECIES}.lnc_RNA.spliced.renamed.fasta"),
    ("lncRNA", "decoy", f"{SPECIES}.lnc_RNA.decoy.spliced.renamed.fasta"),
    ("mRNA", "real", f"{SPECIES}.mRNA.spliced.renamed.fasta"),
    ("mRNA", "decoy", f"{SPECIES}.mRNA.decoy.spliced.renamed.fasta"),
]


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def make_subset(src_fasta, dest_fasta, n, rng):
    """Randomly sample up to `n` records from src_fasta and write them to dest_fasta.

    Returns the number of records written.
    """
    records = list(SeqIO.parse(str(src_fasta), "fasta"))
    if len(records) > n:
        records = rng.sample(records, n)
    SeqIO.write(records, str(dest_fasta), "fasta")
    return len(records)


def tables_to_long(mfe_tables, rna_class, source):
    """Flatten {seq_id: DataFrame(rows=window_size, cols=center)} into a tidy long frame.

    Columns: seq_id, rna_class, source, window_size, center, mfe.
    """
    rows = []
    for seq_id, table in mfe_tables.items():
        # table.index = window sizes; table.columns = center positions
        long = (
            table.reset_index()
            .melt(id_vars="index", var_name="center", value_name="mfe")
            .rename(columns={"index": "window_size"})
            .dropna(subset=["mfe"])
        )
        long["seq_id"] = seq_id
        long["rna_class"] = rna_class
        long["source"] = source
        rows.append(long)
    if not rows:
        return pd.DataFrame(
            columns=["seq_id", "rna_class", "source", "window_size", "center", "mfe"]
        )
    return pd.concat(rows, ignore_index=True)


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
def plot_distribution(df, value, title, out_png):
    """Box+strip distribution of `value` by rna_class x source, faceted by window_size."""
    df = df.copy()
    df["group"] = df["rna_class"] + " / " + df["source"]
    g = sns.catplot(
        data=df,
        x="group",
        y=value,
        col="window_size",
        kind="box",
        order=sorted(df["group"].unique()),
        height=4,
        aspect=0.9,
        showfliers=False,
        color="#9ecae1",
    )
    g.set_xticklabels(rotation=30, ha="right")
    g.set_axis_labels("", value)
    g.figure.suptitle(title, y=1.02)
    g.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(g.figure)


def plot_positional_profile(df, out_png):
    """MFE vs relative window-center position, faceted by window_size.

    Each individual sequence is drawn as a thin, semi-transparent line; the per-group
    mean is overlaid as a thick line.
    """
    df = df.copy()
    df["group"] = df["rna_class"] + " / " + df["source"]
    # Normalise center to [0, 1] within each sequence so transcripts of different
    # lengths can be placed on a common axis.
    df["rel_center"] = df.groupby(["seq_id", "window_size"])["center"].transform(
        lambda c: (c - c.min()) / (c.max() - c.min()) if c.max() > c.min() else 0.0
    )
    df["rel_bin"] = (df["rel_center"] * 20).round() / 20  # 21 bins along 0..1

    # Per-sequence profile (one value per seq x window x bin) and group-mean profile.
    per_seq = (
        df.groupby(["group", "window_size", "seq_id", "rel_bin"], as_index=False)["mfe"]
        .mean()
    )
    prof = (
        df.groupby(["group", "window_size", "rel_bin"], as_index=False)["mfe"].mean()
    )

    groups = sorted(df["group"].unique())
    window_sizes = sorted(df["window_size"].unique())
    palette = dict(zip(groups, sns.color_palette("tab10", len(groups))))

    fig, axes = plt.subplots(
        1, len(window_sizes), figsize=(5.5 * len(window_sizes), 4.5), sharey=False
    )
    if len(window_sizes) == 1:
        axes = [axes]

    for ax, w in zip(axes, window_sizes):
        # thin per-sequence lines
        for (grp, seq_id), sub in per_seq[per_seq["window_size"] == w].groupby(
            ["group", "seq_id"]
        ):
            sub = sub.sort_values("rel_bin")
            ax.plot(
                sub["rel_bin"], sub["mfe"], color=palette[grp],
                lw=0.4, alpha=0.15, zorder=1,
            )
        # thick group-mean lines
        for grp, sub in prof[prof["window_size"] == w].groupby("group"):
            sub = sub.sort_values("rel_bin")
            ax.plot(
                sub["rel_bin"], sub["mfe"], color=palette[grp],
                lw=2.5, alpha=0.95, zorder=3, label=grp,
            )
        ax.set_title(f"window_size = {w}")
        ax.set_xlabel("relative position along transcript (0=5', 1=3')")
        ax.set_ylabel("MFE")

    handles = [
        plt.Line2D([], [], color=palette[g], lw=2.5, label=g) for g in groups
    ]
    axes[-1].legend(handles=handles, title="group", loc="best", fontsize=8)
    fig.suptitle("Positional MFE profile (thin = per sequence, thick = group mean)", y=1.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    for d in (SUBSET_DIR, RESULTS_DIR, PLOTS_DIR):
        d.mkdir(parents=True, exist_ok=True)

    rng = random.Random(SEED)
    all_long = []

    for rna_class, source, fname in INPUTS:
        src = FASTA_DIR / fname
        if not src.exists():
            raise FileNotFoundError(f"Missing input FASTA: {src}")

        subset = SUBSET_DIR / f"{rna_class}_{source}.{N_PER_CLASS}.fasta"
        n = make_subset(src, subset, N_PER_CLASS, rng)
        print(f"[{rna_class}/{source}] sampled {n} sequences -> {subset.name}")

        mfe_tables = helpers.extract_windows_with_mfe_centered(
            str(subset), window_sizes=WINDOW_SIZES, step=STEP
        )
        long = tables_to_long(mfe_tables, rna_class, source)
        all_long.append(long)

    df = pd.concat(all_long, ignore_index=True)
    df["window_size"] = df["window_size"].astype(int)
    df["center"] = df["center"].astype(int)
    df["mfe_per_nt"] = df["mfe"] / df["window_size"]

    # ----- write CSVs ------------------------------------------------------- #
    long_csv = RESULTS_DIR / "mfe_long.csv"
    df.to_csv(long_csv, index=False)
    print(f"Wrote {long_csv} ({len(df):,} rows)")

    summary = (
        df.groupby(["rna_class", "source", "window_size"])
        .agg(
            n_windows=("mfe", "size"),
            n_seqs=("seq_id", "nunique"),
            mfe_mean=("mfe", "mean"),
            mfe_median=("mfe", "median"),
            mfe_std=("mfe", "std"),
            mfe_per_nt_mean=("mfe_per_nt", "mean"),
            mfe_per_nt_median=("mfe_per_nt", "median"),
        )
        .reset_index()
    )
    summary_csv = RESULTS_DIR / "mfe_summary.csv"
    summary.to_csv(summary_csv, index=False)
    print(f"Wrote {summary_csv}")

    # ----- plots ------------------------------------------------------------ #
    sns.set_theme(style="whitegrid")
    plot_distribution(
        df, "mfe", "MFE distribution by class x source", PLOTS_DIR / "mfe_distribution.png"
    )
    plot_distribution(
        df,
        "mfe_per_nt",
        "Length-normalised MFE (MFE / window size)",
        PLOTS_DIR / "mfe_per_nt_distribution.png",
    )
    plot_positional_profile(df, PLOTS_DIR / "positional_profile.png")
    print(f"Wrote 3 plots to {PLOTS_DIR}")

    print("\nSummary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
