#!/usr/bin/env python3
"""
Extract a curated set of metrics from `samtools stats` output and emit
a one-row MultiQC custom-content TSV.

Columns (in order):
  raw_total_sequences, reads_mapped, reads_mapped_percent,
  non_primary_alignments, error_rate, average_quality, average_length

`samtools stats` writes SN summary lines like:
    SN  raw total sequences:    1234
    SN  reads mapped:           5678
    SN  reads mapped %:         98.50
    SN  non-primary alignments: 123
    SN  error rate:             1.23e-3
    SN  average quality:        38.5
    SN  average length:         100
"""
import argparse
import sys
from typing import Dict, List, Tuple


# (output column, raw SN key in samtools stats output)
METRICS: List[Tuple[str, str]] = [
    ("raw_total_sequences",     "raw total sequences"),
    ("reads_mapped",            "reads mapped"),
    ("reads_mapped_percent",    "reads mapped %"),
    ("non_primary_alignments",  "non-primary alignments"),
    ("error_rate",              "error rate"),
    ("average_quality",         "average quality"),
    ("average_length",          "average length"),
]


def parse_stats(path: str) -> Dict[str, str]:
    """Parse SN lines from samtools stats output → {raw_key: value}."""
    out: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("SN\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            key = parts[1].rstrip(":").strip()
            # value may have trailing inline comment after a `#`
            val = parts[2].split("#", 1)[0].strip()
            out[key] = val
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",  required=True, metavar="STATS",
                   help="`samtools stats` output file")
    p.add_argument("--name",   required=True, metavar="SAMPLE",
                   help="Row sample identifier")
    p.add_argument("--output", required=True, metavar="TSV",
                   help="Output one-row TSV")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    data = parse_stats(args.input)

    # `samtools stats` doesn't emit `reads mapped %` directly — derive it.
    try:
        total  = float(data["raw total sequences"])
        mapped = float(data["reads mapped"])
        if total > 0:
            data["reads mapped %"] = f"{(mapped / total) * 100:.2f}"
    except (KeyError, ValueError):
        pass

    row = [data.get(raw, "NA") for _col, raw in METRICS]
    cols = [col for col, _raw in METRICS]

    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write("Sample\t" + "\t".join(cols) + "\n")
        fh.write(args.name + "\t" + "\t".join(row) + "\n")

    n_resolved = sum(1 for v in row if v != "NA")
    print(
        f"SAMTOOLS_TO_MQC: {args.name} → {n_resolved}/{len(METRICS)} metrics resolved",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
