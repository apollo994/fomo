#!/usr/bin/env python3
"""
Drop transcripts (and their lineage) from a GFF3 by a list of transcript IDs.

Used to remove TD2 coding-potential lncRNA from a GFF3. Handles both the flat
projected GFF3 (gene `gene-<tid>` → transcript `<tid>` → `<tid>.exonN`) and the
hierarchical Ensembl-style input GFF3 (gene `gene:<gid>` → transcript
`transcript:<tid>` → exon/CDS). Matching is on a NORMALISED bare id: a leading
`transcript:`, `gene:`, `gene-`, `rna-`, `mrna-`, `cds-` prefix is stripped from
both the GFF feature IDs and the drop-list entries before comparison.

Lineage-aware: a transcript is dropped if its normalised ID is in the drop list;
its child exon/CDS/UTR features (linked via Parent) are dropped with it; a gene
is dropped only if ALL of its transcript children are dropped (so shared genes
with surviving isoforms are kept).

KEEP mode (--keep-fasta + --id-map, plans/20_longest_isoform.md) is the inverse,
used by FILTER_ANNOTATION to build a source's filtered GFF3 FROM the FASTA that is
actually projected: keep exactly the transcripts whose record is in the FASTA,
drop every other one (same lineage rules). The id map (RENAME_FASTA_HEADERS)
translates each renamed FASTA header back to its GFF3 transcript ID verbatim, so
matching is EXACT — normalise() is not used. Fails loudly (exit 1) rather than
letting the GFF3 and the FASTA describe different sets:
  * a renamed header appears twice in the map (two source ids renamed alike);
  * a FASTA record is not in the map;
  * a mapped id is not a transcript of the GFF3.
On success the transcripts of the output == the records of the FASTA.

Usage: filter_gff_by_id.py --gff in.gff3 --drop-ids ids.txt --output out.gff3
       filter_gff_by_id.py --gff in.gff3 --keep-fasta kept.fa --id-map map.tsv --output out.gff3
"""
import argparse
import re
import sys
from typing import Dict, Optional, Set

GENE_TYPES = {"gene", "ncrna_gene", "pseudogene"}
TX_TYPES = {
    "transcript", "mrna", "lnc_rna", "lncrna", "ncrna", "mirna", "trna",
    "rrna", "snrna", "snorna", "pseudogenic_transcript", "primary_transcript",
}
PREFIX_RE = re.compile(r"^(transcript:|gene:|gene-|rna-|mrna-|cds-|cds:)")


def normalise(raw_id: str) -> str:
    return PREFIX_RE.sub("", raw_id.strip())


def attr(col9: str, key: str) -> str:
    """Return the value of an attribute key from a GFF3 column 9, or ''."""
    for field in col9.strip().split(";"):
        field = field.strip()
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return ""


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gff", required=True, help="input GFF3")
    p.add_argument("--drop-ids", help="newline-separated transcript IDs to drop (drop mode)")
    p.add_argument("--keep-fasta", help="FASTA whose records are the transcripts to keep (keep mode)")
    p.add_argument("--id-map", help="TSV: renamed FASTA header -> GFF3 transcript ID (keep mode)")
    p.add_argument("--output", required=True, help="output filtered GFF3")
    args = p.parse_args()
    keep_mode = args.keep_fasta is not None or args.id_map is not None
    if keep_mode and args.drop_ids is not None:
        p.error("--drop-ids is mutually exclusive with --keep-fasta/--id-map")
    if keep_mode and (args.keep_fasta is None or args.id_map is None):
        p.error("keep mode needs both --keep-fasta and --id-map")
    if not keep_mode and args.drop_ids is None:
        p.error("give either --drop-ids or --keep-fasta + --id-map")
    return args


def fail(msg: str) -> None:
    print(f"filter_gff_by_id: ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def read_keep_ids(keep_fasta: str, id_map: str) -> Set[str]:
    """GFF3 transcript IDs (verbatim) of the records in keep_fasta, via id_map."""
    mapping: Dict[str, str] = {}
    with open(id_map, encoding="utf-8") as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                fail(f"{id_map}:{n}: expected 2 tab-separated columns, got {len(parts)}")
            renamed, original = parts
            if renamed in mapping:
                fail(f"renamed header {renamed!r} maps to both {mapping[renamed]!r} and "
                     f"{original!r} — two source transcript ids collapse to one after "
                     f"RENAME_FASTA_HEADERS")
            mapping[renamed] = original
    keep: Set[str] = set()
    with open(keep_fasta, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name not in mapping:
                    fail(f"FASTA record {name!r} is not in the id map {id_map}")
                keep.add(mapping[name])
    return keep


def transcript_ids(raw_lines) -> Set[str]:
    ids: Set[str] = set()
    for line in raw_lines:
        if line.startswith("#") or not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) >= 9 and cols[2].lower() in TX_TYPES:
            tid = attr(cols[8], "ID")
            if tid:
                ids.add(tid)
    return ids


def main() -> int:
    args = parse_args()

    with open(args.gff, encoding="utf-8") as fh:
        raw_lines = fh.readlines()

    # Drop mode matches on normalised ids; keep mode turns the keep set into an
    # exact (verbatim) drop set, so both share the lineage logic below.
    drop: Set[str] = set()
    exact_drop: Optional[Set[str]] = None
    if args.drop_ids is not None:
        with open(args.drop_ids, encoding="utf-8") as fh:
            drop = {normalise(line) for line in fh if line.strip()}
    else:
        keep = read_keep_ids(args.keep_fasta, args.id_map)
        all_tx = transcript_ids(raw_lines)
        missing = sorted(keep - all_tx)
        if missing:
            fail(f"{len(missing)} FASTA record(s) have no transcript in {args.gff}, "
                 f"e.g. {missing[0]!r}")
        exact_drop = all_tx - keep

    # Pass 1: identify dropped transcripts and tally per-gene transcript fate.
    dropped_tx: Set[str] = set()          # raw transcript IDs to drop
    gene_total: Dict[str, int] = {}        # raw gene ID -> n transcript children
    gene_dropped: Dict[str, int] = {}      # raw gene ID -> n dropped transcript children

    for line in raw_lines:
        if line.startswith("#") or not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        ftype = cols[2].lower()
        if ftype in TX_TYPES:
            tid = attr(cols[8], "ID")
            parents = attr(cols[8], "Parent")
            if not tid:
                is_dropped = False
            elif exact_drop is not None:
                is_dropped = tid in exact_drop
            else:
                is_dropped = normalise(tid) in drop
            if is_dropped and tid:
                dropped_tx.add(tid)
            for parent in parents.split(",") if parents else []:
                parent = parent.strip()
                if not parent:
                    continue
                gene_total[parent] = gene_total.get(parent, 0) + 1
                if is_dropped:
                    gene_dropped[parent] = gene_dropped.get(parent, 0) + 1

    # Genes whose every transcript child was dropped.
    dropped_genes = {
        g for g, total in gene_total.items()
        if total > 0 and gene_dropped.get(g, 0) == total
    }

    # Pass 2: emit, skipping dropped transcripts, their children, and dead genes.
    n_dropped = 0
    with open(args.output, "w", encoding="utf-8") as out:
        for line in raw_lines:
            if line.startswith("#") or not line.strip():
                out.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                out.write(line)
                continue
            ftype = cols[2].lower()
            fid = attr(cols[8], "ID")
            parents = [p.strip() for p in attr(cols[8], "Parent").split(",") if p.strip()]

            if ftype in GENE_TYPES:
                if fid in dropped_genes:
                    n_dropped += 1
                    continue
            elif ftype in TX_TYPES:
                if fid in dropped_tx:
                    n_dropped += 1
                    continue
            else:
                # child feature (exon/CDS/UTR): drop if any parent is a dropped tx
                if any(p in dropped_tx for p in parents):
                    n_dropped += 1
                    continue
            out.write(line)

    print(f"filter_gff_by_id: dropped {len(dropped_tx)} transcripts "
          f"({len(dropped_genes)} genes), {n_dropped} feature lines removed", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
