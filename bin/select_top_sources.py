#!/usr/bin/env python3
"""
Select the top-N source species by lncRNA transcript-level F1 from gffcompare
stats, for building a "top-N consensus" annotation.

For each per-source lncRNA gffcompare `.stats` file
(<target>.from_<source>.lnc_RNA.gffcompare.stats), read the Transcript-level
Sensitivity/Precision, compute F1 = 2·Sn·Pr/(Sn+Pr) (0 when Sn+Pr == 0), rank
sources descending, and write the top-N source IDs as a one-column CSV.

Input : one or more <...>.from_<source>.lnc_RNA[.decoy].gffcompare.stats files
Output: --output CSV with header `source` and up to N rows.
"""
import argparse
import sys
from pathlib import Path
from typing import Dict

# Nextflow only puts bin/ on PATH, not on Python's import path — but it stages the
# whole directory, so a __file__-relative insert finds the sibling module. The
# identity grammar, AGGREGATE_IDS and the stats parser live there because this
# script and gffcompare_accuracy_mqc.py must agree on all three.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fomo_stats import AGGREGATE_IDS, f1, parse_identity, transcript_sn_pr  # noqa: E402


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
        # Rank pool is real, individual-source lncRNA only — exclude decoys and
        # the consensus pseudo-sources (the all-source 'combined' and any 'top3').
        if ident.feature_type != "lnc_RNA" or ident.decoy or ident.source in AGGREGATE_IDS:
            continue
        sn_pr = transcript_sn_pr(path)
        if sn_pr is None:
            print(f"WARNING: no Transcript level in '{path}', scoring 0", file=sys.stderr)
            scores[ident.source] = 0.0
            continue
        scores[ident.source] = f1(*sn_pr)

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
