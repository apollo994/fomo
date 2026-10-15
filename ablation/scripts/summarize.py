#!/usr/bin/env python3
"""Turn the per-target ablation tables (ablation/tables/*.tsv) into the Markdown tables quoted
in ablation/ABLATION_REPORT.md. Medians are over the 23 benchmarked targets unless stated.

Usage: summarize.py TABLES_DIR > TABLES_DIR/summary_tables.md
"""
import collections
import csv
import statistics as st
import sys

D = sys.argv[1]


def rd(name):
    return list(csv.DictReader(open(f'{D}/{name}'), delimiter='\t'))


def med(rows, col):
    v = [float(r[col]) for r in rows if r[col] not in ('NA', '')]
    return st.median(v) if v else float('nan')


def tot(rows, col):
    return sum(int(r[col]) for r in rows if r[col] not in ('NA', ''))


def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join('---:' if i else ':---' for i in range(len(header))) + '|']
    for r in rows:
        out.append('| ' + ' | '.join(f'{x:,.1f}' if isinstance(x, float) else (f'{x:,}' if isinstance(x, int) else str(x))
                                     for x in r) + ' |')
    return '\n'.join(out)


def pct(a, b):
    return 100.0 * a / b if b else float('nan')


print('# Ablation summary tables\n')

# ---------------------------------------------------------------- A1 threshold sweep
sw = rd('A1_sweep.tsv')
for ft in ('lnc_RNA', 'mRNA'):
    rows = []
    for k in range(1, 13):
        R = [r for r in sw if r['feature_type'] == ft and int(r['k']) == k]
        rows.append([k, tot(R, 'n_pred'), med(R, 'sn_strict'), med(R, 'pr_strict'), med(R, 'f1_strict'),
                     pct(tot(R, 'tp_lenient'), tot(R, 'n_pred')), tot(R, 'n_decoy'), tot(R, 'novel')])
    print(f'## A1 support threshold — {ft}\n')
    print(table(['k (min species)', 'supported chains (sum)', 'Sn strict (median)', 'Pr strict (median)',
                 'F1 strict (median)', 'Pr lenient (pooled)', 'decoy chains (sum)', 'curated after ref filter (sum)'], rows))
    print()

# ---------------------------------------------------------------- A2 reference filter
R = [r for r in sw if r['feature_type'] == 'lnc_RNA' and r['k'] == '2']
print('## A2 reference filter (lncRNA, k = 2, sums over targets)\n')
print(table(['supported chains', 'kept (novel)', 'lost to reference', 'lost: = annotated lncRNA isoform',
             'lost: = annotated mRNA isoform', 'lost: touches mRNA exon'],
            [[tot(R, 'n_pred'), tot(R, 'novel'), tot(R, 'lost_to_ref'), tot(R, 'lost_eq_annot_lnc'),
              tot(R, 'lost_eq_annot_mrna'), tot(R, 'lost_hits_mrna')]]))
print()

# ---------------------------------------------------------------- A3 calibration
cal = rd('A3_calib.tsv')
bins = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 6), (7, 9), (10, 14), (15, 99)]
for ft in ('lnc_RNA', 'mRNA'):
    rows = []
    for lo, hi in bins:
        R = [r for r in cal if r['feature_type'] == ft and lo <= int(r['support']) <= hi]
        if not R:
            continue
        c = tot(R, 'chains')
        rows.append([f'{lo}' if lo == hi else f'{lo}–{hi if hi < 99 else "+"}', c, pct(tot(R, 'in_strict'), c),
                     pct(tot(R, 'in_lenient'), c), tot(R, 'decoy_chains')])
    print(f'## A3 calibration by exact support — {ft}\n')
    print(table(['species supporting the chain', 'chains', '% = strict reference', '% = any annotated isoform',
                 'decoy chains'], rows))
    print()

# ---------------------------------------------------------------- A4 donors
dn = rd('A4_donors.tsv')
by = collections.defaultdict(dict)
for r in dn:
    if r['strategy'] != 'random':
        by[r['target']][(r['strategy'], int(r['n_sources']), int(r['k']))] = r
cong = [t for t in by if any(k[0] == 'congeners_only' and k[2] == 2 for k in by[t])]
other = [t for t in by if t not in cong]


def pick(ts, key):
    out = []
    for t in ts:
        d = by[t]
        if callable(key):
            r = key(d)
        else:
            r = d.get(key)
        if r:
            out.append(r)
    return out


for label, ts in ((f'targets with >= 2 congeneric donors (n = {len(cong)})', cong),
                  (f'targets with < 2 congeneric donors (n = {len(other)})', other)):
    specs = [('all 22 donors, k = 1 (union)', ('nearest', 22, 1)),
             ('all 22 donors, k = 2 (pipeline default)', ('nearest', 22, 2)),
             ('nearest donor only, k = 1', ('nearest', 1, 1))]
    specs += [(f'nearest {n}, k = 2', ('nearest', n, 2)) for n in (2, 3, 4, 6, 8, 12)]
    specs += [('congeners only, k = 2', lambda d: next((v for k, v in d.items() if k[0] == 'congeners_only' and k[2] == 2), None)),
              ('no congeners, k = 2', lambda d: next((v for k, v in d.items() if k[0] == 'no_congeners' and k[2] == 2), None)),
              ('farthest 12, k = 2', ('farthest', 12, 2))]
    rows = []
    for name, key in specs:
        R = pick(ts, key)
        if R:
            rows.append([name, len(R), med(R, 'n_pred'), med(R, 'sn_strict'), med(R, 'pr_strict'), med(R, 'f1_strict'),
                         med(R, 'f1_lenient')])
    print(f'## A4 donor sets — {label}\n')
    print(table(['donor set', 'targets', 'chains (median)', 'Sn strict', 'Pr strict', 'F1 strict', 'F1 lenient'], rows))
    print()

rows = []
for n in (1, 2, 3, 4, 6, 8, 12, 16):
    for k in (1, 2):
        R = [r for r in dn if r['strategy'] == 'random' and int(r['n_sources']) == n and int(r['k']) == k]
        if R:
            rows.append([n, k, len(R), med(R, 'sn_strict'), med(R, 'pr_strict'), med(R, 'f1_strict')])
print('## A4 random donor subsets (5 draws per target)\n')
print(table(['donors', 'k', 'target × draw', 'Sn strict', 'Pr strict', 'F1 strict'], rows))
print()

# ---------------------------------------------------------------- A5 divergence
ps = rd('A5_persource.tsv')
edges = [(0, .02), (.02, .04), (.04, .06), (.06, .10), (.10, 1)]
rows = []
for lo, hi in edges:
    row = [f'{lo:.2f}–{hi:.2f}' if hi < 1 else f'≥ {lo:.2f}']
    for ft in ('lnc_RNA', 'mRNA'):
        R = [r for r in ps if r['feature_type'] == ft and r['median_de_mRNA'] != 'NA'
             and lo <= float(r['median_de_mRNA']) < hi]
        f1 = [2 * float(r['sn']) * float(r['pr']) / (float(r['sn']) + float(r['pr'])) if float(r['sn']) + float(r['pr']) else 0
              for r in R]
        row += [len(R), st.median(f1) if f1 else float('nan')]
    rows.append(row[:1] + [row[1], row[2], row[4]])
print('## A5 single-donor accuracy vs divergence (median de of the donor\'s projected mRNA)\n')
print(table(['median de', '(target, donor) pairs', 'lncRNA F1 (median)', 'mRNA F1 (median)'], rows))
print()

# ---------------------------------------------------------------- B configs
bc = rd('B_configs.tsv')
cost = collections.Counter()
for r in rd('B_cost.tsv'):
    cost[r['config']] += float(r['cpu_h'])
cov = {r['config']: r for r in rd('B_coding_overlap.tsv')}
base = {(r['target'], r['k']): r for r in bc if r['config'] == 'B0_baseline'}
cfgs = list(dict.fromkeys(r['config'] for r in bc))
rows = []
for c in cfgs:
    R2 = [r for r in bc if r['config'] == c and r['k'] == '2']
    R1 = [r for r in bc if r['config'] == c and r['k'] == '1']
    B2 = [base[(r['target'], '2')] for r in R2]
    dF1 = [float(r['f1_strict']) - float(b['f1_strict']) for r, b in zip(R2, B2)]
    rows.append([c, tot(R2, 'real_in'), tot(R2, 'own_ref_size'), tot(R2, 'proj_td2_coding'), tot(R2, 'raw_models'),
                 tot(R2, 'raw_single_exon'), tot(R2, 'curated_novel_tx'),
                 med(R2, 'sn_strict'), med(R2, 'pr_strict'), med(R2, 'f1_strict'), f'{st.median(dF1):+.2f}',
                 med(R2, 'f1_lenient'), med(R1, 'f1_strict'), tot(R2, 'eq_annot_mrna'),
                 cov[c]['pct'] if c in cov else '—', tot(R2, 'curated_decoy_tx'), tot(R1, 'n_decoy'),
                 round(cost[c], 2)])
print('## B pipeline configurations (lncRNA, fixed yardsticks, k = 2 unless stated)\n')
print(table(['config', 'sequences projected', 'own reference (sum)', 'projected TD2 coding calls', 'raw models',
             'of which single-exon', 'curated transcripts',
             'Sn strict', 'Pr strict', 'F1 strict', 'ΔF1 vs B0 (median paired)', 'F1 lenient', 'F1 strict, k = 1',
             'chains = annotated mRNA', '% chains on mRNA exon', 'curated decoys', 'decoy chains k = 1',
             'CPU-h re-executed'], rows))
print()

# ---------------------------------------------------------------- appendix: per target
dist = {(r['target'], r['rank']): r for r in rd('A4_donors.distance.tsv')}
rows = []
for t in sorted(by):
    d = by[t]
    nd = dist[(t, '1')]
    g = lambda key: float(d[key]['f1_strict'])
    rows.append([t.replace('_', ' '), nd['source'].replace('_', ' '), f"{float(nd['median_de']):.3f}", g(('nearest', 22, 1)), g(('nearest', 22, 2)),
                 g(('nearest', 1, 1)), g(('nearest', 4, 2)), int(t in cong)])
print('## Appendix — per target (lncRNA F1 strict)\n')
print(table(['target', 'nearest donor (by median de)', 'its median de', 'all 22, k = 1', 'all 22, k = 2 (default)',
             'nearest 1, k = 1', 'nearest 4, k = 2', '>= 2 congeners'], rows))
