#!/usr/bin/env python3
"""Ablation analyses for FOMO — the numbers behind ablation/ABLATION_REPORT.md.

Every analysis scores intron CHAINS (multi-exon, stranded, exact intron coordinates) — the
unit curation works on and the unit gffcompare's transcript level matches on under -M (a
multi-exon '=' match is an identical intron chain). Two yardsticks per target:

  strict   the pipeline's own benchmarking reference: <target>.<ft>.filtered.gff3 (spliced,
           longest isoform per gene, lncRNA TD2-noncoding) of a FIXED run, so that configs
           which change the reference (longest-isoform, input TD2) are scored on one ruler
  lenient  every stranded multi-exon transcript of that type in the samplesheet GFF3 (all
           isoforms, no TD2) — a chain matching ANY annotated isoform counts as correct

Inputs are the chain audit tables written by chains_for_run.sh (curate_models.py
--min-species 1 --chains-tsv, with and without --ref), so any support threshold and any donor
subset is evaluated without re-running curation.

Subcommands
  sweep      support threshold k (curation's --min-species) x reference filter, per target
  calib      precision as a function of a chain's exact support
  donors     donor-set ablation: nearest / farthest / random n donors, congeners removed
  persource  per-(target, source) accuracy vs divergence (pipeline gffcompare + median de)
  configs    Part B: score each pipeline configuration on the fixed yardsticks
"""
import argparse
import collections
import csv
import gzip
import json
import math
import os
import random
import re
import statistics
import sys

ID_RE = re.compile(r'(?:^|;)ID=([^;]+)')
PAR_RE = re.compile(r'(?:^|;)Parent=([^;]+)')
SRC_RE = re.compile(r'(?:^|;)source=([^;]+)')
DE_RE = re.compile(r'(?:^|;)de=([^;]+)')


def opener(p):
    return gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p)


def read_sheet(path):
    rows = list(csv.DictReader(open(path)))
    return {r['species']: r for r in rows}


# --------------------------------------------------------------------------- chains
def gff_chains(path, types):
    """{tid: (key, type)} for stranded multi-exon transcripts of column-3 type in `types`;
    key = (seqid, strand, ((exon_end, next_exon_start), ...)) — curate_models.py's chain."""
    tx, ex = {}, collections.defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line[0] == '#':
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            if f[2] == 'exon':
                m = PAR_RE.search(f[8])
                if m:
                    for p in m.group(1).split(','):
                        ex[p].append((f[0], f[6], int(f[3]), int(f[4])))
            elif f[2] in types:
                m = ID_RE.search(f[8])
                if m:
                    tx[m.group(1)] = (f[0], f[6], f[2])
    out = {}
    for tid, (c, s, t) in tx.items():
        e = ex.get(tid)
        if not e or s not in '+-' or len(e) < 2 or any(x[0] != c or x[1] != s for x in e):
            continue
        e = sorted((a, b) for _, _, a, b in e)
        out[tid] = ((c, s, tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))), t)
    return out


def read_chain_table(path):
    """chains.tsv.gz -> {key: (species frozenset, n_species_post, status, ref_types)}."""
    out = {}
    with gzip.open(path, 'rt') as fh:
        for d in csv.DictReader(fh, delimiter='\t'):
            chain = tuple(tuple(int(x) for x in iv.split('-')) for iv in d['introns'].split(','))
            out[(d['chr'], d['strand'], chain)] = (frozenset(d['species_pre'].split(',')),
                                                   int(d['n_species_post']), d['status'], d['ref_types'])
    return out


def prf(pred, truth):
    tp = len(pred & truth)
    sn = 100.0 * tp / len(truth) if truth else float('nan')
    pr = 100.0 * tp / len(pred) if pred else float('nan')
    f1 = 2 * sn * pr / (sn + pr) if pred and truth and (sn + pr) > 0 else 0.0
    return tp, sn, pr, f1


class Yard:
    """Per-target yardsticks, cached."""

    def __init__(self, fixed_results, sheet):
        self.res, self.sheet, self.cache, self.ann = fixed_results, sheet, {}, {}

    def get(self, target, ft):
        k = (target, ft)
        if k not in self.cache:
            # one target at a time: the yardsticks of 23 annotations do not fit a small job
            self.cache = {kk: v for kk, v in self.cache.items() if kk[0] == target}
            strict = gff_chains(f'{self.res}/sources/{target}/annotation/{target}.{ft}.filtered.gff3', {ft})
            if target not in self.ann:      # the samplesheet GFF3 is big: parse it once per target
                self.ann = {target: gff_chains(self.sheet[target]['gff3'], {'lnc_RNA', 'mRNA'})}
            ann = self.ann[target]
            self.cache[k] = dict(
                strict={v[0] for v in strict.values()},
                lenient={v[0] for v in ann.values() if v[1] == ft},
                mrna={v[0] for v in ann.values() if v[1] == 'mRNA'},
                lnc={v[0] for v in ann.values() if v[1] == 'lnc_RNA'})
        return self.cache[k]


def seqkit_inputs(results, any_target):
    """{(source, gtype): n sequences sent to projection} from a target's MultiQC SeqKit table."""
    p = f'{results}/targets/{any_target}/multiqc/multiqc_report_data/multiqc_seqkit_stats_table.txt'
    out = {}
    for d in csv.DictReader(open(p), delimiter='\t'):
        sp, gt = d['Sample'].split('.', 1)
        out[(sp, gt)] = int(float(d['num_seqs']))
    return out


def chains_path(chains_dir, target, track, mode):
    return f'{chains_dir}/{target}/{target}.{track}.{mode}.chains.tsv.gz'


def targets_of(sheet):
    return sorted(s for s, r in sheet.items() if r['role'] == 'both')


def fmt(x):
    if isinstance(x, float):
        return 'NA' if math.isnan(x) else f'{x:.3f}'
    return str(x)


def write_tsv(path, header, rows):
    with open(path, 'w') as o:
        o.write('\t'.join(header) + '\n')
        for r in rows:
            o.write('\t'.join(fmt(r[h]) for h in header) + '\n')
    print(f'wrote {path} ({len(rows)} rows)', file=sys.stderr)


def score(pred, y, decoy_pred=None, n_real_in=None, n_decoy_in=None):
    r = {}
    for name in ('strict', 'lenient'):
        tp, sn, pr, f1 = prf(pred, y[name])
        r.update({f'tp_{name}': tp, f'sn_{name}': sn, f'pr_{name}': pr, f'f1_{name}': f1})
    r['n_pred'] = len(pred)
    r['eq_annot_mrna'] = len(pred & y['mrna'])          # predicted lncRNA chain == an annotated mRNA isoform
    if decoy_pred is not None:
        r['n_decoy'] = len(decoy_pred)
        r['est_fdr'] = (100.0 * (len(decoy_pred) / n_decoy_in) / (len(pred) / n_real_in)
                        if pred and n_decoy_in and n_real_in else float('nan'))
    return r


# --------------------------------------------------------------------------- sweep
def cmd_sweep(a):
    sheet = read_sheet(a.sheet)
    yard = Yard(a.fixed_results, sheet)
    targets = targets_of(sheet)
    inputs = seqkit_inputs(a.results, targets[0])
    rows = []
    for t in targets:
        for ft in a.feature_types:
            y = yard.get(t, ft)
            real = read_chain_table(chains_path(a.chains, t, ft, 'noref'))
            ref = read_chain_table(chains_path(a.chains, t, ft, 'ref'))
            dec = read_chain_table(chains_path(a.chains, t, f'{ft}.decoy', 'noref'))
            srcs = {s for s in targets if s != t}
            n_real_in = sum(inputs.get((s, ft), 0) for s in srcs)
            n_decoy_in = sum(inputs.get((s, f'{ft}.decoy'), 0) for s in srcs)
            for k in range(1, a.kmax + 1):
                pred = {key for key, v in real.items() if len(v[0]) >= k}
                dpred = {key for key, v in dec.items() if len(v[0]) >= k}
                r = dict(target=t, feature_type=ft, k=k, n_ref_strict=len(y['strict']), n_ref_lenient=len(y['lenient']))
                r.update(score(pred, y, dpred, n_real_in, n_decoy_in))
                # with the reference filter, as the pipeline runs: novel = support counted on
                # models with no exon on a same-strand reference exon
                novel = {key for key, v in ref.items() if v[1] >= k}
                lost = {key for key, v in ref.items() if len(v[0]) >= k and v[1] < k}
                r['novel'] = len(novel)
                r['lost_to_ref'] = len(lost)
                r['lost_eq_annot_lnc'] = len(lost & y['lnc'])
                r['lost_eq_annot_mrna'] = len(lost & y['mrna'])
                r['lost_hits_mrna'] = sum('mRNA' in ref[key][3].split(',') for key in lost)
                r['novel_eq_annot'] = len(novel & (y['lnc'] | y['mrna']))
                rows.append(r)
    hdr = ['target', 'feature_type', 'k', 'n_ref_strict', 'n_ref_lenient', 'n_pred', 'tp_strict', 'sn_strict',
           'pr_strict', 'f1_strict', 'tp_lenient', 'sn_lenient', 'pr_lenient', 'f1_lenient', 'eq_annot_mrna',
           'n_decoy', 'est_fdr', 'novel', 'lost_to_ref', 'lost_eq_annot_lnc', 'lost_eq_annot_mrna',
           'lost_hits_mrna', 'novel_eq_annot']
    write_tsv(a.out, hdr, rows)


# --------------------------------------------------------------------------- calib
def cmd_calib(a):
    sheet = read_sheet(a.sheet)
    yard = Yard(a.fixed_results, sheet)
    targets = targets_of(sheet)
    rows = []
    for ft in a.feature_types:
        tally = collections.defaultdict(lambda: [0, 0, 0, 0, 0])   # n -> chains, strict, lenient, mrna, decoy
        for t in targets:
            y = yard.get(t, ft)
            for key, v in read_chain_table(chains_path(a.chains, t, ft, 'noref')).items():
                c = tally[len(v[0])]
                c[0] += 1
                c[1] += key in y['strict']
                c[2] += key in y['lenient']
                c[3] += key in y['mrna']
            for key, v in read_chain_table(chains_path(a.chains, t, f'{ft}.decoy', 'noref')).items():
                tally[len(v[0])][4] += 1
        for n in sorted(tally):
            c = tally[n]
            rows.append(dict(feature_type=ft, support=n, chains=c[0], in_strict=c[1], in_lenient=c[2], eq_annot_mrna=c[3],
                             decoy_chains=c[4], pct_strict=100.0 * c[1] / c[0] if c[0] else float('nan'),
                             pct_lenient=100.0 * c[2] / c[0] if c[0] else float('nan')))
    write_tsv(a.out, ['feature_type', 'support', 'chains', 'in_strict', 'pct_strict', 'in_lenient', 'pct_lenient',
                      'eq_annot_mrna', 'decoy_chains'], rows)


# --------------------------------------------------------------------------- divergence
def median_de(results, target, track):
    """{source: median de} of the target's projected models of one track (truth-free distance)."""
    vals = collections.defaultdict(list)
    with opener(f'{results}/targets/{target}/annotation/{target}.allModels.{track}.raw.gff3') as fh:
        for line in fh:
            if '\ttranscript\t' not in line:
                continue
            f = line.split('\t')
            s, d = SRC_RE.search(f[8]), DE_RE.search(f[8])
            if s and d and d.group(1) != '.':
                vals[s.group(1)].append(float(d.group(1)))
    return {s: statistics.median(v) for s, v in vals.items()}


def genus(sp):
    return sp.split('_')[0]


# --------------------------------------------------------------------------- donors
def cmd_donors(a):
    sheet = read_sheet(a.sheet)
    yard = Yard(a.fixed_results, sheet)
    targets = targets_of(sheet)
    inputs = seqkit_inputs(a.results, targets[0])
    rng = random.Random(a.seed)
    rows, dist_rows = [], []
    ft = a.feature_type
    for t in targets:
        y = yard.get(t, ft)
        real = read_chain_table(chains_path(a.chains, t, ft, 'noref'))
        dec = read_chain_table(chains_path(a.chains, t, f'{ft}.decoy', 'noref'))
        de = median_de(a.results, t, a.distance_track)
        srcs = sorted((s for s in targets if s != t), key=lambda s: (de.get(s, float('inf')), s))
        for i, s in enumerate(srcs, 1):
            dist_rows.append(dict(target=t, source=s, rank=i, median_de=de.get(s, float('nan')),
                                  same_genus=int(genus(s) == genus(t))))
        subsets = []
        for n in a.sizes:
            if n > len(srcs):
                continue
            subsets.append(('nearest', n, 0, srcs[:n]))
            subsets.append(('farthest', n, 0, srcs[-n:]))
            if n < len(srcs):
                for rep in range(1, a.reps + 1):
                    subsets.append(('random', n, rep, rng.sample(srcs, n)))
        cong = [s for s in srcs if genus(s) == genus(t)]
        if cong:
            other = [s for s in srcs if genus(s) != genus(t)]
            subsets.append(('no_congeners', len(other), 0, other))
            subsets.append(('congeners_only', len(cong), 0, cong))
        for strategy, n, rep, sub in subsets:
            S = set(sub)
            n_real_in = sum(inputs.get((s, ft), 0) for s in S)
            n_decoy_in = sum(inputs.get((s, f'{ft}.decoy'), 0) for s in S)
            for k in (1, 2, 3):
                if k > len(S):
                    continue
                pred = {key for key, v in real.items() if len(v[0] & S) >= k}
                dpred = {key for key, v in dec.items() if len(v[0] & S) >= k}
                r = dict(target=t, strategy=strategy, n_sources=n, rep=rep, k=k,
                         n_congeners=sum(genus(s) == genus(t) for s in S),
                         mean_rank=statistics.mean(srcs.index(s) + 1 for s in S))
                r.update(score(pred, y, dpred, n_real_in, n_decoy_in))
                rows.append(r)
    write_tsv(a.out, ['target', 'strategy', 'n_sources', 'rep', 'k', 'n_congeners', 'mean_rank', 'n_pred',
                      'tp_strict', 'sn_strict', 'pr_strict', 'f1_strict', 'tp_lenient', 'sn_lenient', 'pr_lenient',
                      'f1_lenient', 'eq_annot_mrna', 'n_decoy', 'est_fdr'], rows)
    write_tsv(a.out.replace('.tsv', '.distance.tsv'), ['target', 'source', 'rank', 'median_de', 'same_genus'], dist_rows)


# --------------------------------------------------------------------------- persource
def gffcompare_tx(path):
    out = {}
    for line in open(path):
        if 'level:' in line:
            toks = line.replace('|', '').replace('Intron chain', 'Intron_chain').split()
            try:
                out[toks[0]] = (float(toks[2]), float(toks[3]))
            except (ValueError, IndexError):
                pass
    return out


def cmd_persource(a):
    sheet = read_sheet(a.sheet)
    targets = targets_of(sheet)
    rows = []
    for t in targets:
        de = {tr: median_de(a.results, t, tr) for tr in a.distance_tracks}
        for s in targets:
            if s == t:
                continue
            for ft in a.feature_types:
                p = f'{a.results}/targets/{t}/gffcompare/{t}.from_{s}.{ft}.gffcompare.stats'
                if not os.path.exists(p):
                    continue
                acc = gffcompare_tx(p).get('Transcript', (float('nan'), float('nan')))
                r = dict(target=t, source=s, feature_type=ft, same_genus=int(genus(s) == genus(t)),
                         sn=acc[0], pr=acc[1])
                for tr in a.distance_tracks:
                    r[f'median_de_{tr}'] = de[tr].get(s, float('nan'))
                rows.append(r)
    write_tsv(a.out, ['target', 'source', 'feature_type', 'same_genus', 'sn', 'pr'] +
              [f'median_de_{tr}' for tr in a.distance_tracks], rows)


# --------------------------------------------------------------------------- configs
def count_models(path):
    """(transcripts, single-exon transcripts) in a fomo projected GFF3."""
    n_ex = collections.Counter()
    n_tx = 0
    with opener(path) as fh:
        for line in fh:
            f = line.split('\t')
            if len(f) < 9:
                continue
            if f[2] == 'transcript':
                n_tx += 1
            elif f[2] == 'exon':
                n_ex[PAR_RE.search(f[8]).group(1)] += 1
    return n_tx, sum(1 for v in n_ex.values() if v == 1)


def cmd_configs(a):
    sheet = read_sheet(a.sheet)
    yard = Yard(a.fixed_results, sheet)
    targets = targets_of(sheet)
    rows = []
    for spec in a.runs:
        name, res, chains = spec.split(':')
        inputs = seqkit_inputs(res, targets[0])
        for t in targets:
            # the projected-stage TD2 row is target-scoped: only this target's report has it
            p = f'{res}/targets/{t}/multiqc/multiqc_report_data/multiqc_td2_coding_table.txt'
            td2 = {d['Sample']: d for d in csv.DictReader(open(p), delimiter='\t')}
            ft = 'lnc_RNA'
            y = yard.get(t, ft)
            real = read_chain_table(chains_path(chains, t, ft, 'noref'))
            dec = read_chain_table(chains_path(chains, t, f'{ft}.decoy', 'noref'))
            srcs = {s for s in targets if s != t}
            n_real_in = sum(inputs.get((s, ft), 0) for s in srcs)
            n_decoy_in = sum(inputs.get((s, f'{ft}.decoy'), 0) for s in srcs)
            n_tx, n_single = count_models(f'{res}/targets/{t}/annotation/{t}.allModels.{ft}.raw.gff3')
            own = gffcompare_tx(f'{res}/targets/{t}/gffcompare/{t}.from_allModels_raw.{ft}.gffcompare.stats')
            per_src = []
            for s in srcs:
                p = f'{res}/targets/{t}/gffcompare/{t}.from_{s}.{ft}.gffcompare.stats'
                if os.path.exists(p):
                    sn, pr = gffcompare_tx(p).get('Transcript', (0.0, 0.0))
                    per_src.append((2 * sn * pr / (sn + pr)) if sn + pr else 0.0)
            cur = json.load(open(f'{res}/targets/{t}/curated/{t}.curated.{ft}.curation.json'))
            curd = json.load(open(f'{res}/targets/{t}/curated/{t}.curated.{ft}.decoy.curation.json'))
            own_ref = sum(1 for _ in gff_chains(f'{res}/sources/{t}/annotation/{t}.{ft}.filtered.gff3', {ft}))
            proj = td2.get(f'{t}.allModels', {})
            for k in (1, 2):
                pred = {key for key, v in real.items() if len(v[0]) >= k}
                dpred = {key for key, v in dec.items() if len(v[0]) >= k}
                r = dict(config=name, target=t, k=k, real_in=n_real_in, decoy_in=n_decoy_in,
                         raw_models=n_tx, raw_single_exon=n_single,
                         proj_td2_coding=proj.get('n_coding', 'NA'),
                         own_ref_size=own_ref,
                         own_sn_allmodels=own.get('Transcript', (float('nan'),) * 2)[0],
                         own_pr_allmodels=own.get('Transcript', (float('nan'),) * 2)[1],
                         median_source_f1=statistics.median(per_src) if per_src else float('nan'),
                         curated_novel_tx=cur['novel']['transcripts'], curated_novel_genes=cur['novel']['genes'],
                         curated_decoy_tx=curd['novel']['transcripts'])
                r.update(score(pred, y, dpred, n_real_in, n_decoy_in))
                rows.append(r)
    write_tsv(a.out, ['config', 'target', 'k', 'real_in', 'decoy_in', 'raw_models', 'raw_single_exon',
                      'proj_td2_coding', 'own_ref_size', 'own_sn_allmodels', 'own_pr_allmodels', 'median_source_f1',
                      'curated_novel_tx', 'curated_novel_genes', 'curated_decoy_tx', 'n_pred', 'tp_strict',
                      'sn_strict', 'pr_strict', 'f1_strict', 'tp_lenient', 'sn_lenient', 'pr_lenient', 'f1_lenient',
                      'eq_annot_mrna', 'n_decoy', 'est_fdr'], rows)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    def common(sp, chains=True):
        sp.add_argument('--sheet', required=True)
        sp.add_argument('--results', required=True, help='results dir of the run being analysed')
        sp.add_argument('--fixed-results', help='results dir holding the fixed (strict) yardstick; default --results')
        if chains:
            sp.add_argument('--chains', required=True, help='chains_for_run.sh output dir')
        sp.add_argument('--out', required=True)

    s = sub.add_parser('sweep'); common(s)
    s.add_argument('--feature-types', nargs='+', default=['lnc_RNA'])
    s.add_argument('--kmax', type=int, default=12)
    s = sub.add_parser('calib'); common(s)
    s.add_argument('--feature-types', nargs='+', default=['lnc_RNA'])
    s = sub.add_parser('donors'); common(s)
    s.add_argument('--feature-type', default='lnc_RNA')
    s.add_argument('--distance-track', default='mRNA')
    s.add_argument('--sizes', nargs='+', type=int, default=[1, 2, 3, 4, 6, 8, 12, 16, 22])
    s.add_argument('--reps', type=int, default=5)
    s.add_argument('--seed', type=int, default=1)
    s = sub.add_parser('persource'); common(s, chains=False)
    s.add_argument('--feature-types', nargs='+', default=['lnc_RNA', 'mRNA'])
    s.add_argument('--distance-tracks', nargs='+', default=['mRNA', 'lnc_RNA'])
    s = sub.add_parser('configs')
    s.add_argument('--sheet', required=True)
    s.add_argument('--fixed-results', required=True)
    s.add_argument('--runs', nargs='+', required=True, help='NAME:RESULTS_DIR:CHAINS_DIR')
    s.add_argument('--out', required=True)
    a = p.parse_args()
    if getattr(a, 'fixed_results', None) is None:
        a.fixed_results = a.results
    {'sweep': cmd_sweep, 'calib': cmd_calib, 'donors': cmd_donors, 'persource': cmd_persource,
     'configs': cmd_configs}[a.cmd](a)


if __name__ == '__main__':
    main()
