#!/usr/bin/env python3
"""
Extract NM, AS, de values from transcript lines of transferred GFF3 annotations
produced by minimap_transfer, and emit a TSV with source/target/genetype.

Expects layout:
    <root>/<target>_from_<source>/{pcg,lnc,decoy}/aligned.bam.gff3

Genetype mapping: pcg -> mrna, lnc -> lnc, decoy -> decoy.

Usage:
    02_extract_alignment_metrics.py <minimap_models_dir> [-o output.tsv]
"""

import argparse
import re
import sys
from pathlib import Path

GENETYPE_MAP = {"pcg": "mrna", "lnc": "lnc", "decoy": "decoy"}
ATTR_RE = re.compile(r"(?:^|;)(NM|AS|de)=([^;]+)")


def parse_dir_name(name):
    if "_from_" not in name:
        return None, None
    target, source = name.split("_from_", 1)
    return source, target


def extract_metrics(gff3_path):
    rows = []
    with open(gff3_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = dict(ATTR_RE.findall(fields[8]))
            if not {"NM", "AS", "de"}.issubset(attrs):
                continue
            tid_match = re.search(r"(?:^|;)ID=([^;]+)", fields[8])
            tid = tid_match.group(1) if tid_match else ""
            rows.append((tid, attrs["NM"], attrs["AS"], attrs["de"]))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path, help="minimap_models directory")
    ap.add_argument("-o", "--output", type=Path, default=Path("alignment_metrics.tsv"))
    args = ap.parse_args()

    if not args.root.is_dir():
        sys.exit(f"not a directory: {args.root}")

    n_written = 0
    with args.output.open("w") as out:
        out.write("source\ttarget\tgenetype\ttranscript_id\tNM\tAS\tde\n")
        for sample_dir in sorted(p for p in args.root.iterdir() if p.is_dir()):
            source, target = parse_dir_name(sample_dir.name)
            if source is None:
                print(f"skipping (no _from_): {sample_dir.name}", file=sys.stderr)
                continue
            for sub, genetype in GENETYPE_MAP.items():
                gff3 = sample_dir / sub / "aligned.bam.gff3"
                if not gff3.is_file():
                    print(f"missing: {gff3}", file=sys.stderr)
                    continue
                for tid, nm, as_, de in extract_metrics(gff3):
                    out.write(f"{source}\t{target}\t{genetype}\t{tid}\t{nm}\t{as_}\t{de}\n")
                    n_written += 1

    print(f"wrote {n_written} rows -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
