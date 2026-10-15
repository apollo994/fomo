#!/usr/bin/env python3
"""Compute cost of each ablation run: the tasks a configuration actually re-executed
(cached ones are excluded), by process, from its pipeline_info/execution_trace.txt.
MINIMAP2_ALIGN tasks are split by gene type, read from the task's own .command.sh, since the
trace tag is the pseudo-source id 'allModels' for every track.

Usage: trace_cost.py WORKDIR NAME=TRACE [NAME=TRACE ...] > cost.tsv
"""
import collections
import csv
import glob
import re
import sys


def secs(s):
    t = 0.0
    for v, u in re.findall(r'([\d.]+)(ms|s|m|h|d)', s or ''):
        t += float(v) * {'ms': .001, 's': 1, 'm': 60, 'h': 3600, 'd': 86400}[u]
    return t


def track_of(workdir, h):
    for p in glob.glob(f'{workdir}/{h}*/.command.sh'):
        m = re.search(r'allModels\.(lnc_RNA|mRNA)(\.decoy)?\.raw', open(p).read())
        if m:
            return m.group(1) + (m.group(2) or '')
    return '?'


def main():
    workdir = sys.argv[1]
    print('config\tprocess\ttasks\tcpu_h\trealtime_median_s\trealtime_max_s')
    for spec in sys.argv[2:]:
        name, path = spec.split('=', 1)
        agg = collections.defaultdict(list)
        for r in csv.DictReader(open(path), delimiter='\t'):
            if r['status'] != 'COMPLETED':
                continue
            p = r['process'].split(':')[-1]
            if p == 'MINIMAP2_ALIGN':
                p += f'[{track_of(workdir, r["hash"])}]'
            agg[p].append((secs(r['realtime']), int(r['cpus'] or 1)))
        for p, v in sorted(agg.items()):
            rt = sorted(x[0] for x in v)
            cpu_h = sum(x[0] * x[1] for x in v) / 3600
            print(f'{name}\t{p}\t{len(v)}\t{cpu_h:.3f}\t{rt[len(rt) // 2]:.1f}\t{rt[-1]:.1f}')


if __name__ == '__main__':
    main()
