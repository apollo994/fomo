#!/usr/bin/env python3
"""What-if for UNANNOTATED targets: would masking with the projected mRNA track catch the
coding contamination that TD2 misses?

For every supported lncRNA chain (>= K species, reference filter off), its footprint is the
union of the exons of its supporting models. Truth label: one of its models overlaps a
same-strand ANNOTATED mRNA exon of the target (the reference filter's own verdict, read from the
--ref chain table). Prediction: the footprint overlaps a same-strand exon of a PROJECTED mRNA
model (allModels.mRNA.raw.gff3) that is itself supported by >= M species (M = 1: any projected
mRNA). Reports recall / precision of the mask per M, pooled over targets.

Usage: mrna_mask.py SAMPLESHEET RESULTS CHAINS_DIR K M [M ...] > out.tsv
"""
import bisect
import collections
import csv
import gzip
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ablate import PAR_RE, SRC_RE, opener  # noqa: E402


def models(path):
    """{model id: (chr, strand, [exons]), species}"""
    ex, sp = collections.defaultdict(list), {}
    with opener(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            if f[2] == 'exon':
                ex[PAR_RE.search(f[8]).group(1)].append((f[0], f[6], int(f[3]), int(f[4])))
            elif f[2] == 'transcript':
                tid = f[8].split(';', 1)[0][3:]
                sp[tid] = SRC_RE.search(f[8]).group(1)
    return ex, sp


class Index:
    def __init__(self, ivs):
        d = collections.defaultdict(list)
        for c, s, a, b in ivs:
            d[(c, s)].append((a, b))
        self.d = {}
        for k, v in d.items():
            v.sort()
            self.d[k] = ([x[0] for x in v], [x[1] for x in v], max(x[1] - x[0] for x in v))

    def hit(self, c, s, a, b):
        if (c, s) not in self.d:
            return False
        st, en, w = self.d[(c, s)]
        i = bisect.bisect_left(st, a - w)
        while i < len(st) and st[i] <= b:
            if en[i] >= a:
                return True
            i += 1
        return False


def chain_rows(path):
    with gzip.open(path, 'rt') as fh:
        for d in csv.DictReader(fh, delimiter='\t'):
            yield d


def main():
    sheet, res, chains, K = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    Ms = [int(x) for x in sys.argv[5:]]
    targets = [r['species'] for r in csv.DictReader(open(sheet)) if r['role'] == 'both']
    tally = {m: collections.Counter() for m in Ms}
    for t in targets:
        lnc_ex, _ = models(f'{res}/targets/{t}/annotation/{t}.allModels.lnc_RNA.raw.gff3')
        # mRNA chains and their support, from the mRNA no-ref chain table
        m_support = {}
        for d in chain_rows(f'{chains}/{t}/{t}.mRNA.noref.chains.tsv.gz'):
            n = int(d['n_species_pre'])
            for mid in d['models'].split(','):
                m_support[mid] = n
        m_ex, _ = models(f'{res}/targets/{t}/annotation/{t}.allModels.mRNA.raw.gff3')
        idx = {m: Index([e for mid, v in m_ex.items() if m_support.get(mid, 0) >= m for e in v]) for m in Ms}
        for d in chain_rows(f'{chains}/{t}/{t}.lnc_RNA.ref.chains.tsv.gz'):
            if int(d['n_species_pre']) < K:
                continue
            truth = 'mRNA' in d['ref_types'].split(',')
            foot = [e for mid in d['models'].split(',') for e in lnc_ex[mid]]
            for m in Ms:
                pred = any(idx[m].hit(c, s, a, b) for c, s, a, b in foot)
                tally[m][(truth, pred)] += 1
        print(f'done {t}', file=sys.stderr)
    print('min_mrna_species\tsupported_lnc_chains\tcoding_truth\tflagged\ttp\trecall\tprecision\tflagged_not_coding')
    for m in Ms:
        c = tally[m]
        n = sum(c.values())
        truth = c[(True, True)] + c[(True, False)]
        flagged = c[(True, True)] + c[(False, True)]
        print(f'{m}\t{n}\t{truth}\t{flagged}\t{c[(True, True)]}\t{100.0 * c[(True, True)] / truth:.1f}\t'
              f'{100.0 * c[(True, True)] / flagged if flagged else float("nan"):.1f}\t{c[(False, True)]}')


if __name__ == '__main__':
    main()
