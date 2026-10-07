#!/usr/bin/env python3
"""
Split a projected GFF3 by the source species that produced each model.

Since the single-mapping redesign (plans/15) every source's spliced transcripts are
merged into one FASTA and aligned to a target in ONE minimap2 run, so a target has a
single `allModels.raw.gff3` instead of one GFF3 per source. This script recovers the
per-source view from it.

The key is the `source=` attribute that `bin/bam_to_gff.sh` already stamps on every
`gene` and `transcript` row (it lifts it from the `<tid>|<gene_class>|<source>` read
name that `modules/local/rename_fasta_headers.nf` writes). Child features — `exon`,
CDS, UTR — carry no `source=` and inherit it through `Parent`.

  --name-template '<target>.from_{source}.lnc_RNA.projected.gff3'
      One output file per distinct `source=` value. `{source}` is substituted.
      (The subset mode that built the top-N consensus was removed with it, plans/23.)

Records are written in input order, so a coordinate-sorted input yields
coordinate-sorted, seqid-grouped outputs (which `gff-feature-stats` requires).

Duplicate transcript IDs are a hard error. Per-source GFF3s used to be separate
files, so two species sharing a transcript id never collided; merged into one file a
duplicate `ID=` would make `gffread -w` and the gffcompare collapse silently
misbehave. Source transcript ids are NOT guaranteed unique across species —
seen for real, two Drosophila species both used the same generic "lnc_RNA1412"
id — so `bin/bam_to_gff.sh` makes `ID=` globally unique BY CONSTRUCTION
(`<tid>|<source>`, not bare `<tid>`) rather than assuming upstream ids never
collide. This check should therefore never fire now — but it still fails loudly
rather than corrupting the per-source split if that guarantee is ever broken.
"""
import argparse
import os
import sys
from typing import Dict, TextIO

# Feature types that own a transcript-level ID. Mirrors filter_gff_by_id.py's TX_TYPES
# so both scripts agree on what "a transcript" is; in practice bam_to_gff.sh only ever
# writes `transcript`.
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
