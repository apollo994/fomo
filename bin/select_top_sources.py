#!/usr/bin/env python3
"""
Select the top-N source species by lncRNA transcript-level F1 from gffcompare
stats, for building a "top-N consensus" annotation.

For each per-source lncRNA gffcompare `.stats` file
(<target>.from_<source>.lnc_RNA.gffcompare.stats), read the Transcript-level
Sensitivity/Precision, compute F1 = 2·Sn·Pr/(Sn+Pr) (0 when Sn+Pr == 0), rank
sources descending, and write the top-N source IDs as a one-column CSV.

Self-contained (no sibling imports — Nextflow bin/ scripts only share PATH,
not a Python import path).

Input : one or more <...>.from_<source>.lnc_RNA[.decoy].gffcompare.stats files
Output: --output CSV with header `source` and up to N rows.
"""
import argparse
import os
import re
import sys
from typing import Dict, Optional, Tuple

FNAME_RE = re.compile(
    r"^(?P<target>.+?)\.from_(?P<source>.+?)\.(?P<ft>lnc_RNA|mRNA)(?P<decoy>\.decoy)?$"
)

# Pseudo-source ids standing for an aggregate model rather than a real species. They
# use the same `from_<id>` filename convention as per-source stats, so they must be
# skipped here or the top-N would be picked from models built out of the top-N.
# Kept in sync with AGGREGATE_IDS in subworkflows/local/consensus_top.nf, which
# filters the channel before it ever reaches this script — this is the second layer.
AGGREGATE_IDS = ("allModels_raw", "allModels_collapsed", "top3_raw", "top3_collapsed")


def parse_identity(stats_path: str) -> Optional[Tuple[str, str, bool]]:
    """(source, feature_type, decoy) from a gffcompare stats filename."""
    base = os.path.basename(stats_path)
    base = re.sub(r"\.gffcompare\.stats$", "", base)
    base = re.sub(r"\.stats$", "", base)
    base = re.sub(r"\.gffcompare$", "", base)
    m = FNAME_RE.match(base)
    if not m:
        return None
    return m.group("source"), m.group("ft"), bool(m.group("decoy"))


def transcript_sn_pr(stats_path: str) -> Optional[Tuple[float, float]]:
    """Return (sensitivity, precision) from the Transcript-level row, or None."""
    with open(stats_path, encoding="utf-8") as fh:
        for line in fh:
            if "Transcript level:" not in line:
                continue
            toks = line.replace("|", "").split()
            # e.g. ['Transcript', 'level:', '19.1', '15.6']
            if len(toks) < 4:
                return None
            try:
                return float(toks[2]), float(toks[3])
            except ValueError:
                return None
    return None


def f1(sn: float, pr: float) -> float:
    return 0.0 if (sn + pr) == 0 else 2.0 * sn * pr / (sn + pr)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stats", nargs="+", help="per-source lncRNA gffcompare .stats files")
    p.add_argument("--n", type=int, default=3, help="number of top sources to keep (default 3)")
    p.add_argument("--output", required=True, metavar="CSV", help="output CSV of selected source IDs")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    scores: Dict[str, float] = {}
    for path in args.stats:
        ident = parse_identity(path)
        if ident is None:
            print(f"WARNING: cannot parse identity from '{path}', skipping", file=sys.stderr)
            continue
        source, ft, decoy = ident
        # Rank pool is real, individual-source lncRNA only — exclude decoys and
        # the consensus pseudo-sources (the all-source 'combined' and any 'top3').
        if ft != "lnc_RNA" or decoy or source in AGGREGATE_IDS:
            continue
        sn_pr = transcript_sn_pr(path)
        if sn_pr is None:
            print(f"WARNING: no Transcript level in '{path}', scoring 0", file=sys.stderr)
            scores[source] = 0.0
            continue
        scores[source] = f1(*sn_pr)

    if not scores:
        print("ERROR: no per-source lncRNA stats found to rank", file=sys.stderr)
        return 1

    # Rank by F1 desc, tie-break by source name for determinism.
    ranked = sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))
    top = ranked[: args.n]

    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write("source\n")
        for source, _score in top:
            fh.write(f"{source}\n")

    summary = ", ".join(f"{s}={v:.3f}" for s, v in ranked)
    print(f"select_top_sources: ranked lncRNA F1 [{summary}]", file=sys.stderr)
    print(f"select_top_sources: selected top-{len(top)} → "
          f"{', '.join(s for s, _ in top)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
