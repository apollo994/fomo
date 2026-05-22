#!/usr/bin/env python3
"""
07_bam_aligned_metrics.py

Walk minimap_transfer BAM outputs and emit a per-primary-alignment TSV
combining alignment-length metrics (read_length, aligned_length,
aligned_pct) with the alignment tags NM, AS, de pulled directly from the
BAM.

Expects layout (same as scripts/minimap_transfer/02_extract_alignment_metrics.py):

    <root>/<target>_from_<source>/{pcg,lnc,decoy}/aligned.bam

Genetype mapping: pcg -> mrna, lnc -> lnc, decoy -> decoy.

Only primary alignments are considered. Secondary / supplementary records
are ignored. Unmapped reads are kept by default (read_length from the SEQ
field, aligned_length = 0, aligned_pct = 0.00, empty NM/AS/de); pass
--mapped-only to drop them.

Read length comes from infer_read_length() so hard-clipped primaries
(e.g. with `-Y`) still report the original input transcript length.

Output columns:
    source  target  genetype  transcript_id  read_length  aligned_length  aligned_pct  NM  AS  de

Usage:
    ./07_bam_aligned_metrics.py <minimap_models_dir> [-o output.tsv] [--mapped-only]
"""
import argparse
import sys
from pathlib import Path

import pysam

GENETYPE_MAP = {"pcg": "mrna", "lnc": "lnc", "decoy": "decoy"}


def parse_dir_name(name):
    if "_from_" not in name:
        return None, None
    target, source = name.split("_from_", 1)
    return source, target


def get_tag(rec, tag):
    try:
        return rec.get_tag(tag)
    except KeyError:
        return ""


def extract_metrics_from_bam(bam_path, mapped_only=False):
    seen = set()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary:
                continue
            qname = r.query_name
            if qname is None or qname in seen:
                continue
            seen.add(qname)

            if r.is_unmapped:
                if mapped_only:
                    continue
                rl = r.infer_read_length() or r.query_length or 0
                yield qname, rl, 0, 0.0, "", "", ""
                continue

            rl = r.infer_read_length()
            if rl is None:
                rl = r.query_length or 0
            qaln = r.query_alignment_length or 0
            pct = (100.0 * qaln / rl) if rl else 0.0
            de = get_tag(r, "de")
            de_str = f"{de:.4f}" if isinstance(de, float) else (str(de) if de != "" else "")
            yield qname, rl, qaln, pct, get_tag(r, "NM"), get_tag(r, "AS"), de_str


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("root", type=Path, help="minimap_models directory")
    ap.add_argument(
        "-o", "--output", type=Path, default=Path("bam_aligned_metrics.tsv")
    )
    ap.add_argument(
        "--mapped-only",
        action="store_true",
        help="omit unmapped reads from the output",
    )
    args = ap.parse_args()

    if not args.root.is_dir():
        sys.exit(f"not a directory: {args.root}")

    n_written = 0
    with args.output.open("w") as out:
        out.write(
            "source\ttarget\tgenetype\ttranscript_id\t"
            "read_length\taligned_length\taligned_pct\tNM\tAS\tde\n"
        )
        for sample_dir in sorted(p for p in args.root.iterdir() if p.is_dir()):
            source, target = parse_dir_name(sample_dir.name)
            if source is None:
                print(f"skipping (no _from_): {sample_dir.name}", file=sys.stderr)
                continue
            for sub, genetype in GENETYPE_MAP.items():
                bam = sample_dir / sub / "aligned.bam"
                if not bam.is_file():
                    print(f"missing: {bam}", file=sys.stderr)
                    continue
                for tid, rl, qaln, pct, nm, as_, de in extract_metrics_from_bam(
                    bam, mapped_only=args.mapped_only
                ):
                    out.write(
                        f"{source}\t{target}\t{genetype}\t{tid}\t"
                        f"{rl}\t{qaln}\t{pct:.2f}\t{nm}\t{as_}\t{de}\n"
                    )
                    n_written += 1

    print(f"wrote {n_written} rows -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
