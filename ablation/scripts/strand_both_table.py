#!/usr/bin/env python3
"""Curated lncRNA transcripts with the default same-strand reference filter (the reference run's
own curation JSONs) vs --strand both (strand_both.sh output), per target.

Usage: strand_both_table.py SAMPLESHEET RESULTS_DIR STRAND_BOTH_DIR > A2_strand_both.tsv
"""
import csv
import json
import sys

sheet, res, sb_dir = sys.argv[1:4]
print('target\tcurated_same_strand\tcurated_both_strands\tantisense_exonic_genes_same\tgenes_same')
for r in csv.DictReader(open(sheet)):
    if r['role'] != 'both':
        continue
    t = r['species']
    sb = json.load(open(f'{sb_dir}/{t}.lnc_RNA.strand_both.curation.json'))['novel']
    de = json.load(open(f'{res}/targets/{t}/curated/{t}.curated.lnc_RNA.curation.json'))['novel']
    print(f"{t}\t{de['transcripts']}\t{sb['transcripts']}\t{de['genes_location:antisense_exonic']}\t{de['genes']}")
