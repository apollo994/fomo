#!/usr/bin/env python3
"""Fixture tests for bin/curate_models.py (plans/21, phase A). Needs gffcompare and gffread on PATH.

    python3 tests/curate_models/test_curate_models.py        # exits 1 on the first failed check
"""
import gzip
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, '..', '..', 'bin', 'curate_models.py')
TARGET = 'Tgt_1'
FAILS = []


def check(cond, msg):
    print(('ok    ' if cond else 'FAIL  ') + msg)
    if not cond:
        FAILS.append(msg)


# ------------------------------------------------------------------ fixtures
def model(tid, sp, chrom, strand, exons, de=0.1):
    tid = f'{tid}|{sp}'
    s, e = exons[0][0], exons[-1][1]
    rows = [f'{chrom}\tfomo\tgene\t{s}\t{e}\t1\t{strand}\t.\tID=gene-{tid};source={sp};gene_class=lncRNA',
            f'{chrom}\tfomo\ttranscript\t{s}\t{e}\t1\t{strand}\t.\tID={tid};Parent=gene-{tid};source={sp};gene_class=lncRNA;AS=100;de={de}']
    rows += [f'{chrom}\tfomo\texon\t{a}\t{b}\t1\t{strand}\t.\tID={tid}.exon{i};Parent={tid}' for i, (a, b) in enumerate(exons, 1)]
    return rows


M = []
# 1  chain shared by A and B, different ends, equal exonic length -> rep by lowest de (b1)
M += model('a1', 'A', 'chr1', '+', [(100, 200), (300, 400), (500, 600)], de=0.2)
M += model('b1', 'B', 'chr1', '+', [(150, 200), (300, 400), (500, 650)], de=0.1)
# 2  same chain twice in one species -> unsupported
M += model('a2x', 'A', 'chr1', '+', [(1000, 1100), (1200, 1300)])
M += model('a2y', 'A', 'chr1', '+', [(1010, 1100), (1200, 1310)])
# 3  two supported chains sharing one SJ -> one gene, rep = more SJs (X), 4 species in gene
M += model('a3', 'A', 'chr1', '+', [(2000, 2100), (2200, 2300), (2400, 2500)])
M += model('b3', 'B', 'chr1', '+', [(2000, 2100), (2200, 2300), (2400, 2500)])
M += model('c3', 'C', 'chr1', '+', [(2050, 2100), (2200, 2350)])
M += model('d3', 'D', 'chr1', '+', [(2050, 2100), (2200, 2350)])
# 4  two supported chains bridged only by an unsupported model -> two genes
M += model('a4p', 'A', 'chr1', '+', [(3000, 3100), (3200, 3300)])
M += model('b4p', 'B', 'chr1', '+', [(3000, 3100), (3200, 3300)])
M += model('c4q', 'C', 'chr1', '+', [(3500, 3600), (3700, 3800)])
M += model('d4q', 'D', 'chr1', '+', [(3500, 3600), (3700, 3800)])
M += model('e4', 'E', 'chr1', '+', [(3050, 3100), (3200, 3300), (3500, 3600), (3700, 3750)])
# 8  ties: R2 (3 sp) beats R1 (2 sp) at equal SJ count; they share one SJ
M += model('a8r1', 'A', 'chr1', '+', [(5000, 5100), (5200, 5300), (5400, 5500), (5600, 5700)])
M += model('b8r1', 'B', 'chr1', '+', [(5000, 5100), (5200, 5300), (5400, 5500), (5600, 5700)])
M += model('a8r2', 'A', 'chr1', '+', [(5050, 5100), (5200, 5300), (5450, 5500), (5650, 5700)])
M += model('b8r2', 'B', 'chr1', '+', [(5050, 5100), (5200, 5300), (5450, 5500), (5650, 5700)])
M += model('c8r2', 'C', 'chr1', '+', [(5050, 5100), (5200, 5300), (5450, 5500), (5650, 5700)])
# 8b everything equal -> smallest ID
M += model('m_b', 'A', 'chr1', '-', [(6000, 6100), (6200, 6300)])
M += model('m_a', 'B', 'chr1', '-', [(6000, 6100), (6200, 6300)])
# 9  escaped comma in an ID survives as-is into supporting_models
M += model('x%2Cy', 'B', 'chr1', '-', [(9000, 9100), (9200, 9300)])
M += model('xy', 'C', 'chr1', '-', [(9000, 9100), (9200, 9300)])
# 10 single-exon
M += model('a10', 'A', 'chr1', '+', [(8000, 8300)])
M += model('b10', 'B', 'chr1', '+', [(8000, 8300)])
# 14 self + one species -> unsupported unless --keep-self
M += model('t14', TARGET, 'chr1', '+', [(7000, 7100), (7200, 7300)])
M += model('a14', 'A', 'chr1', '+', [(7000, 7100), (7200, 7300)])
# 5  ref tRNA exon removes b5 -> chain loses support (type tRNA)
M += model('a5', 'A', 'chr2', '+', [(100, 200), (300, 400)])
M += model('b5', 'B', 'chr2', '+', [(100, 200), (300, 720)])
# 6  opposite-strand ref exon -> kept, antisense
M += model('a6', 'A', 'chr2', '-', [(1000, 1100), (1200, 1300)])
M += model('b6', 'B', 'chr2', '-', [(1000, 1100), (1200, 1300)])
# 7  inside a ref lncRNA intron -> kept, within ref span
M += model('a7', 'A', 'chr2', '+', [(2000, 2100), (2200, 2300)])
M += model('b7', 'B', 'chr2', '+', [(2000, 2100), (2200, 2300)])
# 16 two ref types on one chain (mRNA + antisense_RNA) -> primary mRNA
M += model('a16', 'A', 'chr2', '+', [(3000, 3100), (3200, 3300)])
M += model('b16', 'B', 'chr2', '+', [(3000, 3100), (3200, 3300)])
# 16 CDS-only reference model, misc_RNA transcript, childless miRNA
M += model('a16c', 'A', 'chr2', '+', [(4000, 4100), (4200, 4300)])
M += model('b16c', 'B', 'chr2', '+', [(4000, 4100), (4200, 4300)])
M += model('a16m', 'A', 'chr2', '+', [(5000, 5100), (5200, 5300)])
M += model('b16m', 'B', 'chr2', '+', [(5000, 5100), (5200, 5300)])
M += model('a16r', 'A', 'chr2', '+', [(6000, 6100), (6200, 6300)])
M += model('b16r', 'B', 'chr2', '+', [(6000, 6100), (6200, 6300)])
# 17 representative = most SJs first: L (3 SJs, 2 species) beats F (1 SJ, 3 species); they share an SJ
M += model('a17l', 'A', 'chr1', '+', [(10000, 10100), (10200, 10300), (10400, 10500), (10600, 10700)])
M += model('b17l', 'B', 'chr1', '+', [(10000, 10100), (10200, 10300), (10400, 10500), (10600, 10700)])
M += model('c17f', 'C', 'chr1', '+', [(10050, 10100), (10200, 10350)])
M += model('d17f', 'D', 'chr1', '+', [(10050, 10100), (10200, 10350)])
M += model('e17f', 'E', 'chr1', '+', [(10050, 10100), (10200, 10350)])
# 18 ref_location classes on chr5 (reference transcripts in REF / REF_GTF)
LOC = {'a18so|A': 'sense_span_overlap', 'a18ai|A': 'antisense_intronic', 'a18dv|A': 'divergent',
       'a18dm|A': 'divergent', 'a18far|A': 'intergenic', 'a18cv|A': 'intergenic',
       'a6|A': 'antisense_exonic', 'a7|A': 'intronic'}
for sp in ('A', 'B'):
    lo = sp.lower()
    M += model(f'{lo}18so', sp, 'chr5', '+', [(1400, 1450), (1600, 1650)])          # straddles ref span 1500-1750 (+)
    M += model(f'{lo}18ai', sp, 'chr5', '-', [(3600, 3650), (3700, 3750)])          # inside ref (+) intron 3550-3800
    M += model(f'{lo}18dv', sp, 'chr5', '+', [(5000, 5100), (5200, 5300)])          # ref (-) ends 4800: 5' ends 200 bp apart
    M += model(f'{lo}18far', sp, 'chr5', '+', [(8000, 8100), (8200, 8300)])        # ref (-) ends 6990: 1010 bp > 1000
    M += model(f'{lo}18cv', sp, 'chr5', '+', [(13000, 13100), (13200, 13300)])      # ref (-) downstream: tail-to-tail
    M += model(f'{lo}18dm', sp, 'chr5', '-', [(20000, 20100), (20200, 20300)])      # ref (+) starts 20500: 200 bp
# 15 a supported chain on chr3 (in the reference) to test placement in the merged file
M += model('a15', 'A', 'chr3', '+', [(100, 200), (300, 400)])
M += model('b15', 'B', 'chr3', '+', [(100, 200), (300, 400)])

REF = """##gff-version 3
chr2\tRefSeq\tregion\t1\t100000\t.\t+\t.\tID=chr2:1..100000
chr2\tRefSeq\tgene\t700\t750\t.\t+\t.\tID=gene-trna;gene_biotype=tRNA
chr2\tRefSeq\ttRNA\t700\t750\t.\t+\t.\tID=rna-trna;Parent=gene-trna
chr2\tRefSeq\texon\t700\t750\t.\t+\t.\tID=exon-trna;Parent=rna-trna
chr2\tRefSeq\tgene\t1000\t1300\t.\t+\t.\tID=gene-m0;gene_biotype=protein_coding
chr2\tRefSeq\tmRNA\t1050\t1080\t.\t+\t.\tID=rna-m0;Parent=gene-m0
chr2\tRefSeq\texon\t1050\t1080\t.\t+\t.\tID=exon-m0;Parent=rna-m0
chr2\tRefSeq\tgene\t1900\t2600\t.\t+\t.\tID=gene-l7;gene_biotype=lncRNA
chr2\tRefSeq\tlnc_RNA\t1900\t2600\t.\t+\t.\tID=rna-l7;Parent=gene-l7
chr2\tRefSeq\texon\t1900\t1950\t.\t+\t.\tID=exon-l7-1;Parent=rna-l7
chr2\tRefSeq\texon\t2500\t2600\t.\t+\t.\tID=exon-l7-2;Parent=rna-l7
chr2\tRefSeq\tgene\t3050\t3060\t.\t+\t.\tID=gene-m1;gene_biotype=protein_coding
chr2\tRefSeq\tmRNA\t3050\t3060\t.\t+\t.\tID=rna-m1;Parent=gene-m1
chr2\tRefSeq\texon\t3050\t3060\t.\t+\t.\tID=exon-m1;Parent=rna-m1
chr2\tRefSeq\tgene\t3250\t3260\t.\t+\t.\tID=gene-as1;gene_biotype=antisense_RNA
chr2\tRefSeq\tantisense_RNA\t3250\t3260\t.\t+\t.\tID=rna-as1;Parent=gene-as1
chr2\tRefSeq\texon\t3250\t3260\t.\t+\t.\tID=exon-as1;Parent=rna-as1
chr2\tRefSeq\tgene\t4050\t4060\t.\t+\t.\tID=gene-cds;gene_biotype=protein_coding
chr2\tRefSeq\tmRNA\t4050\t4060\t.\t+\t.\tID=rna-cds;Parent=gene-cds
chr2\tRefSeq\tCDS\t4050\t4060\t.\t+\t0\tID=cds-cds;Parent=rna-cds
chr2\tRefSeq\tgene\t5050\t5060\t.\t+\t.\tID=gene-misc;gene_biotype=protein_coding
chr2\tRefSeq\ttranscript\t5050\t5060\t.\t+\t.\tID=rna-misc;Parent=gene-misc;gbkey=misc_RNA
chr2\tRefSeq\texon\t5050\t5060\t.\t+\t.\tID=exon-misc;Parent=rna-misc
chr2\tRefSeq\tgene\t6050\t6060\t.\t+\t.\tID=gene-mir;gene_biotype=miRNA
chr2\tRefSeq\tmiRNA\t6050\t6060\t.\t+\t.\tID=rna-mir;Parent=gene-mir
chr3\tRefSeq\tgene\t5000\t6000\t.\t+\t.\tID=gene-c3;gene_biotype=protein_coding
chr3\tRefSeq\tmRNA\t5000\t6000\t.\t+\t.\tID=rna-c3;Parent=gene-c3
chr3\tRefSeq\texon\t5000\t6000\t.\t+\t.\tID=exon-c3;Parent=rna-c3
chr4\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene-c4;gene_biotype=protein_coding
chr4\tRefSeq\tmRNA\t100\t200\t.\t+\t.\tID=rna-c4;Parent=gene-c4
chr4\tRefSeq\texon\t100\t200\t.\t+\t.\tID=exon-c4;Parent=rna-c4
chr5\tRefSeq\tlnc_RNA\t1500\t1750\t.\t+\t.\tID=rna-so
chr5\tRefSeq\texon\t1500\t1550\t.\t+\t.\tID=exon-so-1;Parent=rna-so
chr5\tRefSeq\texon\t1700\t1750\t.\t+\t.\tID=exon-so-2;Parent=rna-so
chr5\tRefSeq\tmRNA\t3500\t3850\t.\t+\t.\tID=rna-ai
chr5\tRefSeq\texon\t3500\t3550\t.\t+\t.\tID=exon-ai-1;Parent=rna-ai
chr5\tRefSeq\texon\t3800\t3850\t.\t+\t.\tID=exon-ai-2;Parent=rna-ai
chr5\tRefSeq\tmRNA\t4500\t4800\t.\t-\t.\tID=rna-dv
chr5\tRefSeq\texon\t4500\t4600\t.\t-\t.\tID=exon-dv-1;Parent=rna-dv
chr5\tRefSeq\texon\t4700\t4800\t.\t-\t.\tID=exon-dv-2;Parent=rna-dv
chr5\tRefSeq\tlnc_RNA\t6500\t6990\t.\t-\t.\tID=rna-far
chr5\tRefSeq\texon\t6500\t6600\t.\t-\t.\tID=exon-far-1;Parent=rna-far
chr5\tRefSeq\texon\t6800\t6990\t.\t-\t.\tID=exon-far-2;Parent=rna-far
chr5\tRefSeq\tmRNA\t13500\t13800\t.\t-\t.\tID=rna-cv
chr5\tRefSeq\texon\t13500\t13600\t.\t-\t.\tID=exon-cv-1;Parent=rna-cv
chr5\tRefSeq\texon\t13700\t13800\t.\t-\t.\tID=exon-cv-2;Parent=rna-cv
chr5\tRefSeq\tmRNA\t20500\t20800\t.\t+\t.\tID=rna-dm
chr5\tRefSeq\texon\t20500\t20600\t.\t+\t.\tID=exon-dm-1;Parent=rna-dm
chr5\tRefSeq\texon\t20700\t20800\t.\t+\t.\tID=exon-dm-2;Parent=rna-dm
"""

REF_GTF = '\n'.join(
    '\t'.join([c, 'ens', 'exon', s, e, '.', st, '.', f'gene_id "{t}"; transcript_id "{t}"; transcript_biotype "{bt}";'])
    for c, s, e, st, t, bt in [
        ('chr2', '700', '750', '+', 'trna', 'tRNA'), ('chr2', '1050', '1080', '+', 'm0', 'protein_coding'),
        ('chr2', '1900', '1950', '+', 'l7', 'lncRNA'), ('chr2', '2500', '2600', '+', 'l7', 'lncRNA'),
        ('chr2', '3050', '3060', '+', 'm1', 'protein_coding'), ('chr2', '3250', '3260', '+', 'as1', 'antisense'),
        ('chr2', '4050', '4060', '+', 'cds', 'protein_coding'), ('chr2', '5050', '5060', '+', 'misc', 'misc_RNA'),
        ('chr2', '6050', '6060', '+', 'mir', 'miRNA'), ('chr3', '5000', '6000', '+', 'c3', 'protein_coding'),
        ('chr5', '1500', '1550', '+', 'so', 'lncRNA'), ('chr5', '1700', '1750', '+', 'so', 'lncRNA'),
        ('chr5', '3500', '3550', '+', 'ai', 'protein_coding'), ('chr5', '3800', '3850', '+', 'ai', 'protein_coding'),
        ('chr5', '4500', '4600', '-', 'dv', 'protein_coding'), ('chr5', '4700', '4800', '-', 'dv', 'protein_coding'),
        ('chr5', '6500', '6600', '-', 'far', 'lncRNA'), ('chr5', '6800', '6990', '-', 'far', 'lncRNA'),
        ('chr5', '13500', '13600', '-', 'cv', 'protein_coding'), ('chr5', '13700', '13800', '-', 'cv', 'protein_coding'),
        ('chr5', '20500', '20600', '+', 'dm', 'protein_coding'), ('chr5', '20700', '20800', '+', 'dm', 'protein_coding')]) + '\n'


def run(tmp, name, gff_rows, extra, ref=None):
    gff = os.path.join(tmp, f'{TARGET}.allModels.lnc_RNA.raw.gff3')
    with open(gff, 'w') as o:
        o.write('##gff-version 3\n' + ''.join(r + '\n' for r in gff_rows))
    prefix = os.path.join(tmp, name, f'{TARGET}.curated.lnc_RNA')
    cmd = [sys.executable, SCRIPT, '--gff', gff, '--prefix', prefix] + extra + (['--ref', ref] if ref else [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr)
    return r, prefix


def read_gff(path):
    genes, tx = {}, {}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line[0] == '#':
                continue
            f = line.rstrip('\n').split('\t')
            a = dict(kv.split('=', 1) for kv in f[8].split(';'))
            if f[2] == 'gene':
                genes[a['ID']] = (f, a)
            elif f[2] != 'exon':
                tx[a['ID']] = (f, a)
    return genes, tx


def main():
    tmp = tempfile.mkdtemp(prefix='curate_test_')
    ref = os.path.join(tmp, 'ref.gff3')
    with open(ref, 'w') as o:
        o.write(REF)
    ref_gtf = os.path.join(tmp, 'ref.gtf')
    with open(ref_gtf, 'w') as o:
        o.write(REF_GTF)

    # ---- main run: GFF3 reference
    r, pre = run(tmp, 'main', M, ['--gene-prefix', 'FOMO_T'], ref)
    check(r.returncode == 0, 'main run exits 0')
    genes, tx = read_gff(pre + '.gff3.gz')
    rep = {a['ID'] for _, a in tx.values()}
    R = json.load(open(pre + '.curation.json'))
    o, d, n = R['observed'], R['discarded_by_reference'], R['novel']
    check('b1|B' in rep and 'a1|A' not in rep, '1  shared chain kept; equal length -> lowest de (b1)')
    check(not any(x.startswith('a2') for x in rep), '2  same-species duplicate chain is not supported')
    g3 = [a for _, a in tx.values() if a['ID'] in ('a3|A', 'b3|B')]
    check(len(g3) == 1 and not any(x in rep for x in ('c3|C', 'd3|D')), '3  chains sharing one SJ -> one gene, rep has more SJs')
    gene3 = genes[g3[0]['Parent']][1] if g3 else {}
    check(gene3.get('n_chains') == '2' and gene3.get('n_species_gene') == '4', '3  gene reports n_chains=2, n_species_gene=4')
    check(g3 and g3[0]['n_supporting_species'] == '2' and g3[0]['supporting_models'] == 'a3|A,b3|B',
          '3  rep lists its own chain support (2 species, a3|A,b3|B)')
    check(len([x for x in rep if x[:3] in ('a4p', 'b4p', 'c4q', 'd4q')]) == 2, '4  unsupported bridge does not merge two genes')
    check(any(x in rep for x in ('a8r2|A', 'b8r2|B', 'c8r2|C')) and not any(x in rep for x in ('a8r1|A', 'b8r1|B')),
          '8  equal SJ count -> more supporting species wins')
    check('m_a|B' in rep and 'm_b|A' not in rep, '8b full tie -> smallest ID')
    t9 = [a for _, a in tx.values() if a['ID'] in ('x%2Cy|B', 'xy|C')]
    check(t9 and t9[0]['supporting_models'] == 'x%2Cy|B,xy|C', '9  escaped comma in ID kept verbatim in supporting_models')
    check(o['single_exon'] == 2, '10 single-exon models counted and skipped')
    check(o['self_excluded'] == 1 and not any(x.startswith(('t14', 'a14')) for x in rep), '14 self projection gives no support')
    check(not any(x.startswith(('a5', 'b5')) for x in rep), '5  ref exon removes b5 -> chain loses support')
    check(d.get('chains_lost_by_type:tRNA') == 1, '5  lost chain reported under tRNA')
    check('a6|A' in rep or 'b6|B' in rep, '6  opposite-strand ref exon does not remove (same-strand rule)')
    check(n.get('genes_location:antisense_exonic') == 1, '6  counted as antisense_exonic')
    check('a7|A' in rep or 'b7|B' in rep, '7  inside a ref intron is kept')
    check(n.get('genes_location:intronic') == 1, '7  counted as intronic')
    check('a17l|A' in rep and not any(x in rep for x in ('b17l|B', 'c17f|C', 'd17f|D', 'e17f|E')),
          '17 chain with more SJs (3 SJs, 2 sp) beats one shared by more species (1 SJ, 3 sp); tie -> smallest ID')
    t17 = [a for _, a in tx.values() if a['ID'] == 'a17l|A']
    check(t17 and t17[0]['n_supporting_species'] == '2' and t17[0]['n_sj'] == '3' and genes[t17[0]['Parent']][1]['n_chains'] == '2'
          and genes[t17[0]['Parent']][1]['n_species_gene'] == '5',
          '17 rep reports 2 species, 3 SJs; gene n_chains=2, n_species_gene=5')
    loc = {a['ID']: a.get('ref_location') for _, a in tx.values()}
    for rid, want in LOC.items():
        check(loc.get(rid) == want, f'18 {rid} ref_location={want} (got {loc.get(rid)})')
    check(all(genes[a['Parent']][1].get('ref_location') == a.get('ref_location') for _, a in tx.values()),
          '18 ref_location identical on gene and transcript rows')
    far = [a for _, a in tx.values() if a['ID'] == 'a18far|A']
    check(far and far[0].get('ref_distance') == '1010' and far[0].get('ref_nearest') == 'lnc_RNA',
          '18 intergenic just beyond --divergent-dist: ref_nearest=lnc_RNA, ref_distance=1010')
    dv = [a for _, a in tx.values() if a['ID'] == 'a18dv|A']
    check(dv and dv[0].get('ref_distance') == '200', '18 divergent: ref_distance=200')
    check(all(a.get('ref_distance') == '0' for _, a in tx.values() if a.get('ref_location') in
              ('antisense_exonic', 'intronic', 'sense_span_overlap', 'antisense_intronic')), '18 overlapping classes have ref_distance=0')
    loc_tot = sum(v for k, v in n.items() if k.startswith('genes_location:'))
    nloc = {k.split(':', 1)[1]: v for k, v in n.items() if k.startswith('genes_location:')}
    check(loc_tot == n['genes'] and len(nloc) == 7 and nloc['no_reference'] == 0 and nloc['divergent'] == 2,
          f'18 all 7 genes_location classes written, sum to genes ({nloc})')
    check(d.get('chains_lost_by_type:mRNA') == 2 and d.get('chains_lost_by_type:antisense_RNA') == 1,
          '16 mRNA + antisense_RNA chain counted once per type (mRNA also from the CDS-only model)')
    check(d.get('chains_lost_by_primary_type:mRNA') == 2 and 'chains_lost_by_primary_type:antisense_RNA' not in d,
          '16 primary type: mRNA wins over antisense_RNA')
    check(d.get('chains_lost_by_type:transcript:misc_RNA') == 1, '16 RefSeq misc_RNA transcript labelled transcript:misc_RNA')
    check(d.get('chains_lost_by_type:miRNA') == 1, '16 childless miRNA counted by its span')
    tot_primary = sum(v for k, v in d.items() if k.startswith('chains_lost_by_primary_type:'))
    check(tot_primary == d['supported_chains_lost'], '16 primary-type counts add up to supported_chains_lost')
    exp_genes = {'b1|B', 'a3|A', 'a4p|A', 'c4q|C', 'a8r2|A', 'm_a|B', 'x%2Cy|B', 'a6|A', 'a7|A', 'a15|A', 'a17l|A'} | set(LOC)
    check(rep == exp_genes, f'expected representatives exactly ({sorted(rep ^ exp_genes)})')
    check(all(a.get('gene_biotype') == 'fomo_lncRNA' for _, a in genes.values()) and
          all(a.get('transcript_biotype') == 'fomo_lncRNA' for _, a in tx.values()), 'biotype fomo_lncRNA on genes and transcripts')

    # merged file
    merged = pre + '.merged.gff3.gz'
    check(os.path.exists(merged), '15 merged file written')
    lines = gzip.open(merged, 'rt').read().splitlines()
    feat = [ln.split('\t') for ln in lines if ln and ln[0] != '#']
    ids = [kv[3:] for f in feat for kv in f[8].split(';') if kv.startswith('ID=')]
    check(len(ids) == len(set(ids)), '15 merged IDs unique')
    order = [f[0] for f in feat]
    blocks = [c for i, c in enumerate(order) if i == 0 or c != order[i - 1]]
    check(len(blocks) == len(set(blocks)), f'15 merged seqids contiguous ({blocks})')
    check(blocks == ['chr2', 'chr3', 'chr4', 'chr5', 'chr1'], '15 curated chr2/chr3/chr5 genes inside their seqid block, chr1 appended')
    ref_lines = [ln for ln in REF.splitlines() if ln]
    check([ln for ln in lines if '\tfomo\t' not in ln and not ln.startswith('# curate_models')] == ref_lines,
          '15 reference lines verbatim and in order')
    r2 = subprocess.run(['gffread', '-E', '-o', os.devnull, '/dev/stdin'], input='\n'.join(lines) + '\n', capture_output=True, text=True)
    check(r2.returncode == 0, '15 gffread accepts merged')

    # 15 a reference without a ##gff-version line: the merged file gets one, then the reference verbatim
    ref_nh = os.path.join(tmp, 'ref_noheader.gff3')
    open(ref_nh, 'w').write(REF.split('\n', 1)[1])
    r, pre_nh = run(tmp, 'noheader', M, ['--gene-prefix', 'FOMO_T'], ref_nh)
    lines_nh = gzip.open(pre_nh + '.merged.gff3.gz', 'rt').read().splitlines()
    check(r.returncode == 0 and lines_nh[0] == '##gff-version 3' and lines_nh[1].startswith('# curate_models')
          and [ln for ln in lines_nh[2:] if '\tfomo\t' not in ln] == ref_lines[1:],
          '15 reference without ##gff-version: header added, reference lines verbatim')

    # 12 no reference
    r, pre2 = run(tmp, 'noref', M, [])
    check(r.returncode == 0 and not os.path.exists(pre2 + '.merged.gff3.gz'), '12 no reference -> no merged file')
    _, tx2 = read_gff(pre2 + '.gff3.gz')
    check({a['ID'] for _, a in tx2.values()} == exp_genes | {'b5|B', 'a16|A', 'a16c|A', 'a16m|A', 'a16r|A'},
          f'12 without reference every supported gene is kept ({len(tx2)})')
    check(all(a.get('ref_location') == 'no_reference' and 'ref_nearest' not in a for _, a in tx2.values()),
          '12 without reference: ref_location=no_reference, no ref_nearest/ref_distance')
    n2 = json.load(open(pre2 + '.curation.json'))['novel']
    check(n2.get('genes_location:no_reference') == n2['genes'] and n2.get('genes_location:intergenic') == 0,
          '12 every location class written; all genes no_reference')
    # 14 keep-self
    r, pre3 = run(tmp, 'keepself', M, ['--keep-self'])
    _, tx3 = read_gff(pre3 + '.gff3.gz')
    check(any(a['ID'] in ('t14|Tgt_1', 'a14|A') for _, a in tx3.values()), '14 --keep-self lets the target support')
    # 13 GTF reference: same representatives as the GFF3 reference
    r, pre4 = run(tmp, 'gtf', M, ['--gene-prefix', 'FOMO_T'], ref_gtf)
    check(r.returncode == 0, '13 GTF reference run exits 0')
    _, tx4 = read_gff(pre4 + '.gff3.gz')
    check({a['ID'] for _, a in tx4.values()} == exp_genes, '13 GTF reference gives the same curated set')
    loc4 = {a['ID']: a.get('ref_location') for _, a in tx4.values()}
    check(all(loc4.get(k) == v for k, v in LOC.items()), '13 GTF reference gives the same ref_location classes')
    # 11 empty input
    r, pre5 = run(tmp, 'empty', [], [], ref)
    check(r.returncode == 0, '11 empty input exits 0')
    g5, _ = read_gff(pre5 + '.gff3.gz')
    m5 = [ln for ln in gzip.open(pre5 + '.merged.gff3.gz', 'rt').read().splitlines() if not ln.startswith('# curate_models')]
    check(not g5 and m5 == ref_lines, '11 empty input -> no genes, merged == reference')
    # decoy / mRNA biotype from --track
    r, pre6 = run(tmp, 'decoy', M[:10], ['--track', 'lnc_RNA.decoy'])
    g6, _ = read_gff(pre6 + '.gff3.gz')
    check(all(a['gene_biotype'] == 'fomo_decoy_lncRNA' for _, a in g6.values()) and g6, 'biotype fomo_decoy_lncRNA from --track')
    r, pre7 = run(tmp, 'mrna', M[:10], ['--track', 'mRNA'])
    g7, _ = read_gff(pre7 + '.gff3.gz')
    check(all(a['gene_biotype'] == 'fomo_mRNA' for _, a in g7.values()) and g7, 'biotype fomo_mRNA from --track')
    # ID clash with the reference
    clash = REF + 'chr4\tRefSeq\tgene\t300\t400\t.\t+\t.\tID=FOMO_T_G000001\n'
    ref_c = os.path.join(tmp, 'ref_clash.gff3')
    open(ref_c, 'w').write(clash)
    r, _ = run(tmp, 'clash', M, ['--gene-prefix', 'FOMO_T'], ref_c)
    check(r.returncode == 1 and 'already used in the reference' in r.stderr, 'curated ID clashing with the reference exits 1')
    # structural errors
    dup = M[:6] + [M[1]]
    r, _ = run(tmp, 'dup', dup, [])
    check(r.returncode == 1 and 'duplicate transcript ID' in r.stderr, 'duplicate transcript ID exits 1')
    # determinism
    r, pre8 = run(tmp, 'main2', M, ['--gene-prefix', 'FOMO_T'], ref)
    same = all(gzip.open(pre + s, 'rb').read() == gzip.open(pre8 + s, 'rb').read() for s in ('.gff3.gz', '.merged.gff3.gz'))
    check(same, 'rerun gives byte-identical GFF3 content')

    print(f'\n{len(FAILS)} failed' if FAILS else '\nall checks passed')
    sys.exit(1 if FAILS else 0)


if __name__ == '__main__':
    main()
