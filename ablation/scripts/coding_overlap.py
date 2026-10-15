#!/usr/bin/env python3
"""Coding contamination seen by an UNANNOTATED target: supported lncRNA chains (>= K species)
whose models overlap a same-strand annotated mRNA exon — what curation would emit for a pure
target with no gff3 — per configuration, read from the --ref chain tables.

Usage: coding_overlap.py SAMPLESHEET K NAME=CHAINS_DIR [...] > out.tsv
"""
import csv
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ablate import read_chain_table, chains_path  # noqa: E402

sheet, K = sys.argv[1], int(sys.argv[2])
targets = [r['species'] for r in csv.DictReader(open(sheet)) if r['role'] == 'both']
print('config\tsupported\thits_mRNA_exon\tpct')
for spec in sys.argv[3:]:
    name, d = spec.split('=', 1)
    sup = hit = 0
    for t in targets:
        for key, v in read_chain_table(chains_path(d, t, 'lnc_RNA', 'ref')).items():
            if len(v[0]) >= K:
                sup += 1
                hit += 'mRNA' in v[3].split(',')
    print(f'{name}\t{sup}\t{hit}\t{100.0 * hit / sup:.2f}')
