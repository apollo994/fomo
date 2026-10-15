#!/usr/bin/env python3
"""Per-task MINIMAP2_ALIGN wall time, keyed by (target, gene type) read from each task's
.command.sh (the trace tag is the pseudo-source 'allModels' for every task).

Usage: align_times.py WORKDIR NAME=TRACE [NAME=TRACE ...] > align_times.tsv
"""
import csv
import glob
import re
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trace_cost import secs  # noqa: E402

print('config\ttarget\tgtype\trealtime_s\tpeak_rss')
wd = sys.argv[1]
for spec in sys.argv[2:]:
    name, path = spec.split('=', 1)
    for r in csv.DictReader(open(path), delimiter='\t'):
        if not r['process'].endswith('MINIMAP2_ALIGN') or r['status'] != 'COMPLETED':
            continue
        for p in glob.glob(f'{wd}/{r["hash"]}*/.command.sh'):
            m = re.search(r'(\S+)\.allModels\.((?:lnc_RNA|mRNA)(?:\.decoy)?)\.raw', open(p).read())
            if m:
                print(f'{name}\t{m.group(1).split("/")[-1]}\t{m.group(2)}\t{secs(r["realtime"]):.1f}\t{r["peak_rss"]}')
                break
