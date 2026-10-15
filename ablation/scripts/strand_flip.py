#!/usr/bin/env python3
"""Are antisense curated lncRNA chains strand flips? For every supported chain (>= K species,
reference filter off), count those sharing at least one intron — same coordinates, OPPOSITE
strand — with an annotated transcript of the target, and those identical to one when flipped.
Same-strand sharing is reported alongside as the reference point.

Usage: strand_flip.py SAMPLESHEET CHAINS_DIR TRACK K > out.tsv
"""
import csv
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ablate import gff_chains, read_chain_table, chains_path  # noqa: E402

sheet, chains, track, K = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
flip = {'+': '-', '-': '+'}
print('target\tsupported\tsame_strand_shared_intron\topp_strand_shared_intron\topp_strand_identical_chain\t'
      'opp_shared_with_mRNA\topp_shared_with_lncRNA')
for r in csv.DictReader(open(sheet)):
    if r['role'] != 'both':
        continue
    t = r['species']
    ann = gff_chains(r['gff3'], {'lnc_RNA', 'mRNA'})
    introns = {}
    chain_set = set()
    for (c, s, ch), typ in ann.values():
        chain_set.add((c, s, ch))
        for iv in ch:
            introns.setdefault((c, s, iv), set()).add(typ)
    sup = [k for k, v in read_chain_table(chains_path(chains, t, track, 'noref')).items() if len(v[0]) >= K]
    same = opp = ident = opp_m = opp_l = 0
    for c, s, ch in sup:
        o = flip[s]
        same += any((c, s, iv) in introns for iv in ch)
        types = set().union(*[introns.get((c, o, iv), set()) for iv in ch])
        opp += bool(types)
        opp_m += 'mRNA' in types
        opp_l += 'lnc_RNA' in types
        ident += (c, o, ch) in chain_set
    print(f'{t}\t{len(sup)}\t{same}\t{opp}\t{ident}\t{opp_m}\t{opp_l}')
