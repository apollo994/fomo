#!/usr/bin/env python3
"""
Split (or subset) a projected GFF3 by the source species that produced each model.

Since the single-mapping redesign (plans/15) every source's spliced transcripts are
merged into one FASTA and aligned to a target in ONE minimap2 run, so a target has a
single `allModels.raw.gff3` instead of one GFF3 per source. This script recovers the
per-source view from it.

The key is the `source=` attribute that `bin/bam_to_gff.sh` already stamps on every
`gene` and `transcript` row (it lifts it from the `<tid>|<gene_class>|<source>` read
name that `modules/local/rename_fasta_headers.nf` writes). Child features — `exon`,
CDS, UTR — carry no `source=` and inherit it through `Parent`.

Two modes:

  split   --name-template '<target>.from_{source}.lnc_RNA.projected.gff3'
          One output file per distinct `source=` value. `{source}` is substituted.

  subset  --keep top_sources.csv --name '<target>.top3.lnc_RNA.raw.gff3'
          One output file holding only the sources listed in the CSV (one column,
          header `source` — the file bin/select_top_sources.py emits).

Records are written in input order, so a coordinate-sorted input yields
coordinate-sorted, seqid-grouped outputs (which `gff-feature-stats` requires).

Duplicate transcript IDs are a hard error. Per-source GFF3s used to be separate
files, so two species sharing a transcript id never collided; merged into one file a
duplicate `ID=` would make `gffread -w` and the gffcompare collapse silently
misbehave. Ensembl ids are species-scoped so this should never fire — but it fails
loudly rather than corrupting the consensus if it ever does.
"""
import argparse
import csv
import os
import sys
from typing import Dict, List, Optional, Set, TextIO

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


def read_keep(path: str) -> Set[str]:
    """Read the one-column `source` CSV emitted by bin/select_top_sources.py."""
    keep: Set[str] = set()
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            source = (row.get("source") or "").strip()
            if source:
                keep.add(source)
    return keep


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff", required=True, help="input GFF3 (an allModels.raw.gff3)")
    p.add_argument("--outdir", required=True, help="directory to write output GFF3s into")
    p.add_argument("--name-template",
                   help="split mode: output filename with a literal {source} placeholder")
    p.add_argument("--keep", help="subset mode: CSV of sources to keep (header `source`)")
    p.add_argument("--name", help="subset mode: output filename")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    subset = args.keep is not None
    if subset:
        if not args.name:
            print("ERROR: --keep requires --name", file=sys.stderr)
            return 2
        keep = read_keep(args.keep)
        if not keep:
            print(f"ERROR: no sources listed in {args.keep}", file=sys.stderr)
            return 1
    elif not args.name_template:
        print("ERROR: one of --name-template (split) or --keep/--name (subset) is required",
              file=sys.stderr)
        return 2
    elif "{source}" not in args.name_template:
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
    subset_handle: Optional[TextIO] = None
    skipped: Dict[str, int] = {}

    def handle_for(source: str) -> Optional[TextIO]:
        """Output handle for a source, or None if this record is not wanted."""
        nonlocal subset_handle
        if subset:
            if source not in keep:
                return None
            if subset_handle is None:
                subset_handle = open(os.path.join(args.outdir, args.name), "w",
                                     encoding="utf-8")
                subset_handle.write(GFF_HEADER)
            return subset_handle
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

                out = handle_for(source)
                if out is None:
                    skipped[source] = skipped.get(source, 0) + 1
                    continue
                out.write(line)
                counts[source] = counts.get(source, 0) + 1
    finally:
        for fh_out in handles.values():
            fh_out.close()
        if subset_handle is not None:
            subset_handle.close()

    if subset:
        missing = sorted(keep - set(counts))
        if missing:
            print(f"WARNING: requested sources absent from {args.gff}: "
                  f"{', '.join(missing)}", file=sys.stderr)
        if subset_handle is None:
            # Every requested source projected nothing. Still emit a valid empty GFF3
            # so downstream gffcompare/gff-feature-stats have something to read.
            with open(os.path.join(args.outdir, args.name), "w", encoding="utf-8") as fh_out:
                fh_out.write(GFF_HEADER)
    elif not handles:
        print(f"WARNING: {args.gff} contained no records — no output files written",
              file=sys.stderr)

    mode = "subset" if subset else "split"
    for source in sorted(set(counts) | set(skipped)):
        kept = counts.get(source, 0)
        dropped = skipped.get(source, 0)
        note = f" (skipped {dropped})" if dropped else ""
        print(f"gff_by_source [{mode}]: {source}\t{kept} records{note}", file=sys.stderr)
    print(f"gff_by_source [{mode}]: {len(source_by_tid)} transcripts across "
          f"{len(set(counts) | set(skipped))} sources", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
