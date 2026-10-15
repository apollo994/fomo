#!/usr/bin/env python3
"""
Filter the transcripts of a GFF3, with their lineage, by transcript ID.

Drop mode (--drop-ids, optional): drop the listed transcripts. Matching is on a normalised bare
id: a leading `transcript:`, `gene:`, `gene-`, `rna-`, `mrna-`, `cds-` or `cds:` is stripped from
both the GFF3 IDs and the list, so the list may come from either naming. Without a list nothing
is dropped (the GFF3 is copied).

Keep mode (--fasta + --id-map): keep exactly the transcripts whose record is in the FASTA. The id
map (RENAME_FASTA_HEADERS) translates each FASTA header back to its GFF3 transcript ID verbatim,
so matching is exact. Exits 1 if a header maps twice, a FASTA record is not in the map, or a
mapped id is not a transcript of the GFF3: on success the output's transcripts are the FASTA's
records.

Lineage, in both modes: a dropped transcript takes its child features (exon, CDS, UTR, ...) with
it, and a gene is dropped only when all of its transcripts are. Every other line is copied.

Usage: filter_gff_by_id.py --gff in.gff3 [--drop-ids ids.txt] --output out.gff3
       filter_gff_by_id.py --gff in.gff3 --fasta kept.fa --id-map map.tsv --output out.gff3
"""

import argparse
import collections
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fomo_io  # noqa: E402

PREFIX_RE = re.compile(r'^(transcript:|gene:|gene-|rna-|mrna-|cds-|cds:)')
log = fomo_io.get_logger('filter_gff_by_id')


def normalise(raw_id: str) -> str:
    return PREFIX_RE.sub('', raw_id.strip())


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--gff', required=True, help='input GFF3')
    p.add_argument('--drop-ids', help='newline-separated transcript IDs to drop (drop mode)')
    p.add_argument('--fasta', help='FASTA whose records are the transcripts to keep (keep mode)')
    p.add_argument('--id-map', help='TSV: renamed FASTA header -> GFF3 transcript ID (keep mode)')
    p.add_argument('--output', required=True, help='output GFF3')
    args = p.parse_args()
    keep_mode = args.fasta is not None or args.id_map is not None
    if keep_mode and args.drop_ids is not None:
        p.error('--drop-ids is mutually exclusive with --fasta/--id-map')
    if keep_mode and (args.fasta is None or args.id_map is None):
        p.error('keep mode needs both --fasta and --id-map')
    return args


def read_id_map(path: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with fomo_io.open_text(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 2:
                fomo_io.die(f'{path}:{n}: expected 2 tab-separated columns, got {len(parts)}')
            renamed, original = parts
            if renamed in mapping:
                fomo_io.die(
                    f'renamed header {renamed!r} maps to both {mapping[renamed]!r} and {original!r} '
                    f'— two source transcript ids collapse to one after RENAME_FASTA_HEADERS'
                )
            mapping[renamed] = original
    return mapping


def read_keep_ids(fasta: str, id_map: str) -> set[str]:
    """GFF3 transcript IDs (verbatim) of the records in the FASTA."""
    mapping = read_id_map(id_map)
    keep: set[str] = set()
    for name in fomo_io.fasta_names(fasta):
        if name not in mapping:
            fomo_io.die(f'FASTA record {name!r} is not in the id map {id_map}')
        keep.add(mapping[name])
    return keep


def read_transcripts(gff: str) -> list[tuple[str, list[str]]]:
    """(ID, parents) of every transcript row, in file order; the ID is '' when the row has none."""
    rows = []
    for _, cols in fomo_io.iter_gff(gff):
        if fomo_io.is_transcript(cols[2]):
            attrs = fomo_io.parse_attrs(cols[8])
            rows.append((attrs.get('ID', ''), fomo_io.split_ids(attrs.get('Parent', ''))))
    return rows


def drop_predicate(args: argparse.Namespace, transcripts: list[tuple[str, list[str]]]):
    """-> is_dropped(transcript ID) for the mode the arguments select."""
    if args.fasta is not None:
        keep = read_keep_ids(args.fasta, args.id_map)
        missing = sorted(keep - {tid for tid, _ in transcripts if tid})
        if missing:
            fomo_io.die(f'{len(missing)} FASTA record(s) have no transcript in {args.gff}, e.g. {missing[0]!r}')
        return lambda tid: tid not in keep
    drop: set[str] = set()
    if args.drop_ids is not None:
        with fomo_io.open_text(args.drop_ids) as fh:
            drop = {normalise(line) for line in fh if line.strip()}
    return lambda tid: normalise(tid) in drop


def select(transcripts: list[tuple[str, list[str]]], is_dropped) -> tuple[set[str], set[str]]:
    """-> (dropped transcript IDs, genes whose every transcript is dropped). A transcript without
    an ID is never dropped, but counts for its genes."""
    dropped_tx: set[str] = set()
    n_total: collections.Counter = collections.Counter()
    n_dropped: collections.Counter = collections.Counter()
    for tid, parents in transcripts:
        dropped = bool(tid) and is_dropped(tid)
        if dropped:
            dropped_tx.add(tid)
        n_total.update(parents)
        if dropped:
            n_dropped.update(parents)
    dropped_genes = {g for g, total in n_total.items() if n_dropped[g] == total}
    return dropped_tx, dropped_genes


def is_dropped_line(cols: list[str], dropped_tx: set[str], dropped_genes: set[str]) -> bool:
    attrs = fomo_io.parse_attrs(cols[8])
    kind = fomo_io.feature_kind(cols[2])
    if kind == 'gene':
        return attrs.get('ID', '') in dropped_genes
    if kind == 'transcript':
        return attrs.get('ID', '') in dropped_tx
    return any(p in dropped_tx for p in fomo_io.split_ids(attrs.get('Parent', '')))


def write_filtered(gff: str, output: str, dropped_tx: set[str], dropped_genes: set[str]) -> int:
    """Copy every line that is not dropped; -> number of feature lines removed."""
    n_removed = 0
    with fomo_io.open_text(gff) as fh, fomo_io.open_text(output, 'wt') as out:
        for line in fh:
            if not line.startswith('#') and line.strip():
                cols = line.rstrip('\n').split('\t')
                if len(cols) >= 9 and is_dropped_line(cols, dropped_tx, dropped_genes):
                    n_removed += 1
                    continue
            out.write(line)
    return n_removed


def main() -> int:
    args = parse_args()
    transcripts = read_transcripts(args.gff)
    dropped_tx, dropped_genes = select(transcripts, drop_predicate(args, transcripts))
    n_removed = write_filtered(args.gff, args.output, dropped_tx, dropped_genes)
    log.info(f'dropped {len(dropped_tx)} transcripts ({len(dropped_genes)} genes), {n_removed} feature lines removed')
    return 0


if __name__ == '__main__':
    fomo_io.run(main)
