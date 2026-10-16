#!/usr/bin/env python3
"""Fraction of a target's ANNOTATED lncRNAs whose exons overlap an opposite-strand exon of
another annotated transcript (any biotype) — the reference rate for curation's
`antisense_exonic` class, to compare the curated genes against.

Usage: antisense_rate.py SAMPLESHEET RESULTS_DIR > antisense.tsv
"""
import bisect
import collections
import csv
import gzip
import re
import sys

PAR = re.compile(r'(?:^|;)Parent=([^;]+)')
ID = re.compile(r'(?:^|;)ID=([^;]+)')


def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def exons(path):
    """{tid: (chr, strand, [(s, e)])} for every exon parent; plus {tid: type}."""
    ex, typ = collections.defaultdict(list), {}
    with opener(path) as fh:
        for line in fh:
            if line[0] == '#':
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            if f[2] == 'exon':
                m = PAR.search(f[8])
                if m:
                    for p in m.group(1).split(','):
                        ex[p].append((f[0], f[6], int(f[3]), int(f[4])))
            else:
                m = ID.search(f[8])
                if m:
                    typ[m.group(1)] = f[2]
    return ex, typ


class Index:
    def __init__(self, ex):
        d = collections.defaultdict(list)
        for tid, v in ex.items():
            for c, s, a, b in v:
                d[(c, s)].append((a, b, tid))
        self.d = {}
        for k, v in d.items():
            v.sort()
            self.d[k] = ([x[0] for x in v], v, max(x[1] - x[0] for x in v))

    def hits(self, c, s, a, b):
        if (c, s) not in self.d:
            return set()
        starts, v, w = self.d[(c, s)]
        i = bisect.bisect_left(starts, a - w)
        out = set()
        while i < len(v) and v[i][0] <= b:
            if v[i][1] >= a:
                out.add(v[i][2])
            i += 1
        return out


def main():
    sheet, res = sys.argv[1], sys.argv[2]
    print('target\tannotated_lncRNA\tantisense_exonic\tpct')
    for r in csv.DictReader(open(sheet)):
        if r['role'] != 'both':
            continue
        t = r['species']
        ex, typ = exons(r['gff3'])
        idx = Index(ex)
        lnc, _ = exons(f'{res}/sources/{t}/annotation/{t}.lnc_RNA.filtered.gff3')
        n = anti = 0
        for tid, v in lnc.items():
            n += 1
            opp = {'+': '-', '-': '+'}.get(v[0][1])
            if opp and any(idx.hits(c, opp, a, b) for c, _, a, b in v):
                anti += 1
        print(f'{t}\t{n}\t{anti}\t{100.0 * anti / n:.1f}')


if __name__ == '__main__':
    main()
