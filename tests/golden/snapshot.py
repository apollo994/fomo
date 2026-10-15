#!/usr/bin/env python3
"""Snapshot the deliverables of a FOMO results directory: one line per published file with a
sha256 of its normalised content (plans/27_refactor.md §5).

Usage: snapshot.py OUTDIR [--info INFO.tsv] > snapshot.tsv

Normalisation (a deliverable is compared on what it means, not on incidental bytes):
  *.bam                      samtools view records, no header (SAMTOOLS env var overrides the binary)
  gffcompare/*.stats         lines not starting with '#' (the header holds the command line)
  *.gff3.gz / *.gz           decompressed bytes
  *.curation.json            canonical JSON without params.gff / params.ref (staging paths)
  run_summary.json           canonical JSON
  everything else            bytes
Excluded: */multiqc/**, pipeline_info/**, *.csi, summary/tables/*_mqc.json (plot specs).
--info writes, for every MultiQC data table, its row count and columns (informational only).
"""

import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


def excluded(rel: Path) -> bool:
    parts = rel.parts
    return (
        'multiqc' in parts
        or parts[0] == 'pipeline_info'
        or rel.suffix == '.csi'
        or (parts[:2] == ('summary', 'tables') and rel.name.endswith('_mqc.json'))
    )


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def normalise(path: Path) -> tuple:
    name = path.name
    if name.endswith('.bam'):
        samtools = os.environ.get('SAMTOOLS', 'samtools')
        out = subprocess.run([samtools, 'view', str(path)], capture_output=True, check=True).stdout
        return 'bam_records', digest(out)
    if path.parent.name == 'gffcompare' and name.endswith('.stats'):
        lines = [ln for ln in path.read_bytes().splitlines(keepends=True) if not ln.startswith(b'#')]
        return 'no_comment_lines', digest(b''.join(lines))
    if name.endswith('.curation.json'):
        doc = json.loads(path.read_text(encoding='utf-8'))
        for key in ('gff', 'ref'):
            doc.get('params', {}).pop(key, None)
        return 'json_no_paths', digest(json.dumps(doc, sort_keys=True).encode())
    if name == 'run_summary.json':
        doc = json.loads(path.read_text(encoding='utf-8'))
        return 'json', digest(json.dumps(doc, sort_keys=True).encode())
    if name.endswith('.gz'):
        with gzip.open(path, 'rb') as fh:
            return 'gunzip', digest(fh.read())
    return 'bytes', digest(path.read_bytes())


def info_rows(outdir: Path):
    tables = (p for p in outdir.rglob('*.txt') if 'multiqc' in p.parts and p.parent.name.endswith('_data'))
    for table in sorted(tables):
        lines = table.read_text(encoding='utf-8', errors='replace').splitlines()
        header = lines[0].split('\t') if lines else []
        yield str(table.relative_to(outdir)), max(len(lines) - 1, 0), ','.join(header)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('outdir', type=Path)
    ap.add_argument('--info', type=Path, help='also write the informational MultiQC table summary here')
    args = ap.parse_args()
    out = sys.stdout
    out.write('path\trule\tsha256\n')
    for path in sorted(p for p in args.outdir.rglob('*') if p.is_file()):
        rel = path.relative_to(args.outdir)
        if excluded(rel):
            continue
        rule, sha = normalise(path)
        out.write(f'{rel}\t{rule}\t{sha}\n')
    if args.info:
        with open(args.info, 'w', encoding='utf-8') as fh:
            fh.write('table\trows\tcolumns\n')
            for rel, rows, cols in info_rows(args.outdir):
                fh.write(f'{rel}\t{rows}\t{cols}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
