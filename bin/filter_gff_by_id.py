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

Usage: filter_gff_by_id.py --gff in.gff3 --drop-ids ids.txt --output out.gff3
"""
import argparse
import re
import sys
from typing import Dict, Set

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
    p.add_argument("--drop-ids", required=True, help="newline-separated transcript IDs to drop")
    p.add_argument("--output", required=True, help="output filtered GFF3")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    with open(args.drop_ids, encoding="utf-8") as fh:
        drop = {normalise(line) for line in fh if line.strip()}

    with open(args.gff, encoding="utf-8") as fh:
        raw_lines = fh.readlines()

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
            is_dropped = normalise(tid) in drop if tid else False
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
