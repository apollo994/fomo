#!/usr/bin/env python3
"""
Split a transferred GFF3 into TP and FP files using gffcompare class codes.

Strict TP definition (per project decision):
    TP = class_code in {'=', 'c', 'k'}
    FP = every other class_code (j, m, n, o, i, y, x, s, u, ...).

Inputs:
    1. <input.gff3>        -- the transferred GFF3 (e.g. benchmark/input/<tool>_<...>.gff3)
    2. <annotated.gtf>     -- gffcompare's <prefix>.annotated.gtf
    3. <out.TP.gff3>       -- where to write TP transcripts (+ their gene/exon/CDS/UTR records)
    4. <out.FP.gff3>       -- where to write FP transcripts

Records handled:
    - 'gene'                                    -> emitted to whichever bucket(s) receive a child transcript
    - 'transcript' | 'mRNA' | 'lnc_RNA' | 'rna' -> classified via class_code
    - all other feature lines (exon, CDS, UTRs) -> follow their Parent transcript

Transcripts not present in the .annotated.gtf are skipped with a warning count.
"""

import re
import sys


TP_CODES = {'=', 'c', 'k'}
TX_TYPES = {'transcript', 'mRNA', 'lnc_RNA', 'rna'}

# matches both `transcript_id "..."` and `class_code "..."` on the same line
CLASS_RE = re.compile(r'transcript_id "([^"]+)".*?class_code "([^"]+)"')


def load_class_codes(annotated_gtf):
    cls = {}
    with open(annotated_gtf) as fh:
        for line in fh:
            if line.startswith('#') or '\ttranscript\t' not in line:
                continue
            m = CLASS_RE.search(line)
            if m:
                cls[m.group(1)] = m.group(2)
    return cls


def parse_attrs(s):
    out = {}
    for part in s.strip().split(';'):
        if '=' in part:
            k, v = part.split('=', 1)
            out[k.strip()] = v.strip()
    return out


def split(gff_in, classes, tp_out, fp_out):
    gene_lines = {}                                   # gene_id -> raw gene line
    gene_emitted = {'TP': set(), 'FP': set()}         # genes already written per bucket
    tx_bucket = {}                                    # tid -> 'TP' | 'FP' (for child lookup)

    out = {'TP': open(tp_out, 'w'), 'FP': open(fp_out, 'w')}
    for f in out.values():
        f.write('##gff-version 3\n')

    n_tp = n_fp = n_skip_tx = 0

    with open(gff_in) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            ftype = fields[2]
            attrs = parse_attrs(fields[8])

            if ftype == 'gene':
                gene_lines[attrs.get('ID', '')] = line
                continue

            if ftype in TX_TYPES:
                tid = attrs.get('ID', '')
                cls = classes.get(tid)
                if cls is None:
                    n_skip_tx += 1
                    continue
                bucket = 'TP' if cls in TP_CODES else 'FP'
                tx_bucket[tid] = bucket
                gid = attrs.get('Parent', '').split(',')[0]
                if gid in gene_lines and gid not in gene_emitted[bucket]:
                    out[bucket].write(gene_lines[gid])
                    gene_emitted[bucket].add(gid)
                out[bucket].write(line)
                if bucket == 'TP':
                    n_tp += 1
                else:
                    n_fp += 1
                continue

            # leaf: exon / CDS / UTR — follow Parent
            parent = attrs.get('Parent', '').split(',')[0]
            bucket = tx_bucket.get(parent)
            if bucket is None:
                continue
            out[bucket].write(line)

    for f in out.values():
        f.close()

    print(f'TP transcripts: {n_tp}', file=sys.stderr)
    print(f'FP transcripts: {n_fp}', file=sys.stderr)
    if n_skip_tx:
        print(f'skipped (no class_code in annotated.gtf): {n_skip_tx}', file=sys.stderr)


def main():
    if len(sys.argv) != 5:
        sys.exit('usage: 05_split_TPFP.py <input.gff3> <annotated.gtf> <out.TP.gff3> <out.FP.gff3>')
    in_gff, annot_gtf, tp_path, fp_path = sys.argv[1:]
    classes = load_class_codes(annot_gtf)
    print(f'loaded {len(classes)} class codes from {annot_gtf}', file=sys.stderr)
    split(in_gff, classes, tp_path, fp_path)


if __name__ == '__main__':
    main()
