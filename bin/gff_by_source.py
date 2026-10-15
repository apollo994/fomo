#!/usr/bin/env python3
"""
Split a projected GFF3 (a target's allModels.raw.gff3, every source in one file) by the source
species that produced each model.

The key is the `source=` attribute that bin/bam_to_gff.sh writes on every gene and transcript row;
child features (exon, CDS, UTR) inherit it through `Parent`.

  --name-template '<target>.from_{source}.lnc_RNA.projected.gff3'
      One output file per distinct `source=` value; `{source}` is substituted.

Records keep their input order, so a coordinate-sorted input gives coordinate-sorted,
seqid-grouped outputs (as gff-feature-stats requires). A duplicate transcript ID is an error:
bam_to_gff.sh makes IDs unique (<tid>|<source>), and a duplicate would corrupt the split and make
gffread and gffcompare misbehave. An input without records writes no file.
"""

import argparse
import collections
import contextlib
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fomo_io  # noqa: E402

GFF_HEADER = '##gff-version 3\n'
log = fomo_io.get_logger('gff_by_source')


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--gff', required=True, help='input GFF3 (an allModels.raw.gff3)')
    p.add_argument('--outdir', required=True, help='directory to write the output GFF3s into')
    p.add_argument('--name-template', required=True, help='output file name with a literal {source} placeholder')
    args = p.parse_args()
    if '{source}' not in args.name_template:
        p.error('--name-template must contain the literal {source}')
    return args


def attributed_rows(gff: str):
    """Yield (source, columns) per feature. Rows carrying `source=` define it; the others take it
    from the first Parent seen so far, which suffices because bam_to_gff.sh writes gene ->
    transcript -> exons in that order."""
    source_by_id: dict[str, str] = {}
    source_by_tid: dict[str, str] = {}
    for lineno, cols in fomo_io.iter_gff(gff):
        attrs = fomo_io.parse_attrs(cols[8])
        fid = attrs.get('ID', '')
        source = attrs.get('source', '') or next(
            (source_by_id[p] for p in fomo_io.split_ids(attrs.get('Parent', '')) if p in source_by_id), ''
        )
        if not source:
            fomo_io.die(
                f'{gff}:{lineno}: cannot attribute record to a source — no source= attribute and no known Parent'
            )
        if fid:
            source_by_id[fid] = source
        if fid and fomo_io.is_transcript(cols[2]):
            if fid in source_by_tid:
                fomo_io.die(
                    f"duplicate transcript ID '{fid}' in {gff} (line {lineno}): claimed by both "
                    f"'{source_by_tid[fid]}' and '{source}'. Transcript IDs must be unique across the merged "
                    f'annotation — gffread and gffcompare silently misbehave otherwise.'
                )
            source_by_tid[fid] = source
        yield source, cols


def split_by_source(gff: str, outdir: str, template: str) -> collections.Counter:
    """Write one GFF3 per source; -> records per source."""
    counts: collections.Counter = collections.Counter()
    os.makedirs(outdir, exist_ok=True)
    with contextlib.ExitStack() as stack:
        handles = {}
        for source, cols in attributed_rows(gff):
            if source not in handles:
                path = os.path.join(outdir, template.replace('{source}', source))
                handles[source] = stack.enter_context(fomo_io.open_text(path, 'wt'))
                handles[source].write(GFF_HEADER)
            handles[source].write('\t'.join(cols) + '\n')
            counts[source] += 1
    return counts


def main() -> int:
    args = parse_args()
    counts = split_by_source(args.gff, args.outdir, args.name_template)
    if not counts:
        log.warning(f'{args.gff} contained no records — no output files written')
    for source in sorted(counts):
        log.info(f'{source}\t{counts[source]} records')
    log.info(f'{sum(counts.values())} records across {len(counts)} sources')
    return 0


if __name__ == '__main__':
    fomo_io.run(main)
