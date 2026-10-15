#!/usr/bin/env python3
"""
Split a projected GFF3 (a target's allModels.raw.gff3, every source in one file) by the
source species that produced each model.

The key is the `source=` attribute that bin/bam_to_gff.sh writes on every `gene` and
`transcript` row; child features (exon, CDS, UTR) inherit it through `Parent`.

  --name-template '<target>.from_{source}.lnc_RNA.projected.gff3'
      One output file per distinct `source=` value; `{source}` is substituted.

Records keep their input order, so a coordinate-sorted input gives coordinate-sorted,
seqid-grouped outputs (as gff-feature-stats requires). A duplicate transcript ID is an
error: bam_to_gff.sh makes IDs unique (<tid>|<source>), and a duplicate would corrupt the
split and make gffread and gffcompare misbehave.
"""
import argparse
import os
import sys
from typing import Dict, TextIO

# Feature types that own a transcript-level ID (the same set as filter_gff_by_id.py).
TX_TYPES = {
    "transcript", "mrna", "lnc_rna", "lncrna", "ncrna", "mirna", "trna",
    "rrna", "snrna", "snorna", "pseudogenic_transcript", "primary_transcript",
}

GFF_HEADER = "##gff-version 3\n"


def attr(col9: str, key: str) -> str:
    """Return the value of an attribute key from a GFF3 column 9, or ''."""
    for field in col9.strip().split(";"):
        field = field.strip()
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return ""


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff", required=True, help="input GFF3 (an allModels.raw.gff3)")
    p.add_argument("--outdir", required=True, help="directory to write output GFF3s into")
    p.add_argument("--name-template", required=True,
                   help="output filename with a literal {source} placeholder")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if "{source}" not in args.name_template:
        print("ERROR: --name-template must contain the literal {source}", file=sys.stderr)
        return 2

    os.makedirs(args.outdir, exist_ok=True)

    # ID -> source, populated as records are read. Rows carrying `source=` seed it;
    # rows without one resolve through Parent. A forward single pass suffices because
    # bam_to_gff.sh emits gene → transcript → exons in that order.
    source_by_id: Dict[str, str] = {}
    # Transcript ID -> the source that first claimed it (for the duplicate guard).
    source_by_tid: Dict[str, str] = {}

    handles: Dict[str, TextIO] = {}
    counts: Dict[str, int] = {}

    def handle_for(source: str) -> TextIO:
        if source not in handles:
            fname = args.name_template.replace("{source}", source)
            handles[source] = open(os.path.join(args.outdir, fname), "w", encoding="utf-8")
            handles[source].write(GFF_HEADER)
        return handles[source]

    try:
        with open(args.gff, encoding="utf-8") as fh:
            for lineno, line in enumerate(fh, 1):
                # Header/comment lines are regenerated per output, not copied: an
                # output must start with its own ##gff-version 3 and nothing else.
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue

                ftype = cols[2].lower()
                fid = attr(cols[8], "ID")
                source = attr(cols[8], "source")

                if not source:
                    # Child feature: inherit from the first resolvable parent.
                    for parent in (p.strip() for p in attr(cols[8], "Parent").split(",")):
                        if parent in source_by_id:
                            source = source_by_id[parent]
                            break

                if not source:
                    print(f"ERROR: {args.gff}:{lineno}: cannot attribute record to a "
                          f"source — no source= attribute and no known Parent",
                          file=sys.stderr)
                    return 1

                if fid:
                    source_by_id[fid] = source

                if ftype in TX_TYPES and fid:
                    previous = source_by_tid.get(fid)
                    if previous is not None:
                        print(f"ERROR: duplicate transcript ID '{fid}' in {args.gff} "
                              f"(line {lineno}): claimed by both '{previous}' and "
                              f"'{source}'. Transcript IDs must be unique across the "
                              f"merged annotation — gffread and gffcompare silently "
                              f"misbehave otherwise.", file=sys.stderr)
                        return 1
                    source_by_tid[fid] = source

                handle_for(source).write(line)
                counts[source] = counts.get(source, 0) + 1
    finally:
        for fh_out in handles.values():
            fh_out.close()

    if not handles:
        print(f"WARNING: {args.gff} contained no records — no output files written",
              file=sys.stderr)

    for source in sorted(counts):
        print(f"gff_by_source: {source}\t{counts[source]} records", file=sys.stderr)
    print(f"gff_by_source: {len(source_by_tid)} transcripts across {len(counts)} sources",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
