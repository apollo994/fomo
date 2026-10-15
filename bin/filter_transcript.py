#!/usr/bin/env python3
"""
Keep the stranded, spliced transcripts of one type, with their parent genes and child features.

Stranded: the transcript's own strand is + or -, every child (exon/CDS/UTR) is on that same
strand, and no parent has an undefined strand. Anything else (trans-spliced models, e.g. RefSeq
plant mitochondrial genes with strand '?') is dropped: gffread -w rejects strand '?', and a
mixed-strand model is neither one cDNA nor one alignment.

Spliced: at least two exons, so the source projections and the benchmarking reference compare
like with like (gffcompare runs with -M, and decoys require two exons too).

--longest keeps one transcript per gene: among the gene's spliced, stranded transcripts of the
type, the one with the greatest summed exon length (what gffread -w extracts and minimap2
aligns). Selection runs after the spliced filter, so a gene whose longest isoform is mono-exonic
keeps its longest spliced one. A gene is the transcript's whole Parent value; a transcript without
Parent competes only with itself; ties go to the smallest transcript ID in byte order.
Rationale: plans/20.

--counts FILE --sample NAME write a one-row funnel table (*_transcript_filter_mqc.tsv):
  n_transcripts  transcripts of the type
  n_stranded     of those, consistently stranded
  n_spliced      of those, with >= 2 exons
  n_selected     written (== n_spliced without --longest)

Usage: filter_transcript.py --gff in.gff3[.gz] --feature-type lnc_RNA [--longest]
                            [--counts FILE --sample NAME] --output out.gff3
"""

import argparse
import collections
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fomo_io  # noqa: E402

STRANDS = ('+', '-')
FUNNEL = ('n_transcripts', 'n_stranded', 'n_spliced', 'n_selected')


@dataclass
class Scan:
    """What the first pass learns about the candidates and their children."""

    parent_of: dict = field(default_factory=dict)  # candidate ID -> its Parent value ('' if none)
    strand: dict = field(default_factory=dict)  # candidate ID -> its own strand
    bad_strand: set = field(default_factory=set)  # IDs of features with an undefined strand
    child_strands: dict = field(default_factory=lambda: collections.defaultdict(set))  # parent -> child strands
    exon_count: collections.Counter = field(default_factory=collections.Counter)
    exon_len: collections.Counter = field(default_factory=collections.Counter)


def scan(gff: str, feature_type: str) -> Scan:
    s = Scan()
    for _, cols in fomo_io.iter_gff(gff):
        attrs = fomo_io.parse_attrs(cols[8])
        fid, parent = attrs.get('ID', ''), attrs.get('Parent', '')
        if fid and cols[6] not in STRANDS:
            s.bad_strand.add(fid)
        if cols[2] == feature_type:
            if fid:
                s.parent_of[fid] = parent
                s.strand[fid] = cols[6]
            continue
        is_exon = cols[2].lower() == 'exon'
        for p in (x for x in parent.split(',') if x):
            s.child_strands[p].add(cols[6] if cols[6] in STRANDS else '?')
            if is_exon:
                s.exon_count[p] += 1
                s.exon_len[p] += int(cols[4]) - int(cols[3]) + 1
    return s


def is_stranded(s: Scan, tid: str) -> bool:
    strand = s.strand[tid]
    if strand not in STRANDS or not s.child_strands.get(tid, set()) <= {strand}:
        return False
    return not any(p in s.bad_strand for p in s.parent_of[tid].split(',') if p)


def select(s: Scan, longest: bool) -> tuple[set[str], dict]:
    """-> (selected transcript IDs, funnel counts)."""
    stranded = [tid for tid in s.parent_of if is_stranded(s, tid)]
    spliced = [tid for tid in stranded if s.exon_count[tid] >= 2]
    if longest:
        best: dict[str, str] = {}
        for tid in spliced:
            gene = s.parent_of[tid] or tid
            other = best.get(gene)
            if other is None or (s.exon_len[tid], other) > (s.exon_len[other], tid):
                best[gene] = tid
        selected = set(best.values())
    else:
        selected = set(spliced)
    counts = dict(zip(FUNNEL, (len(s.parent_of), len(stranded), len(spliced), len(selected)), strict=True))
    return selected, counts


def write_selected(gff: str, output: str, selected: set[str], parents: set[str]) -> None:
    """Write the selected transcripts, their parents and their children, in input order."""
    with fomo_io.open_text(output, 'wt') as out:
        for _, cols in fomo_io.iter_gff(gff):
            attrs = fomo_io.parse_attrs(cols[8])
            fid = attrs.get('ID', '')
            if (fid and (fid in parents or fid in selected)) or any(
                p in selected for p in attrs.get('Parent', '').split(',')
            ):
                out.write('\t'.join(cols) + '\n')


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--gff', required=True, help='input GFF3 (.gz ok)')
    p.add_argument('--feature-type', required=True, help='transcript type to keep, e.g. lnc_RNA')
    p.add_argument('--longest', action='store_true', help='keep the longest spliced transcript per gene')
    p.add_argument('--counts', help='funnel table to write (needs --sample)')
    p.add_argument('--sample', help='row name of the funnel table')
    p.add_argument('--output', required=True, help='output GFF3')
    args = p.parse_args()
    if args.counts and not args.sample:
        p.error('--counts requires --sample')
    return args


def main() -> int:
    args = parse_args()
    s = scan(args.gff, args.feature_type)
    selected, counts = select(s, args.longest)
    parents = {p for tid in selected for p in s.parent_of[tid].split(',') if p}
    write_selected(args.gff, args.output, selected, parents)
    if args.counts:
        with fomo_io.open_text(args.counts, 'wt') as out:
            out.write('\t'.join(('Sample', 'feature') + FUNNEL) + '\n')
            out.write('\t'.join([args.sample, args.feature_type] + [str(counts[k]) for k in FUNNEL]) + '\n')
    return 0


if __name__ == '__main__':
    fomo_io.run(main)
