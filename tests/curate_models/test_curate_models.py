"""pytest suite for bin/curate_models.py: fixture models covering support, ranking, reference filter,
location classes, reference tags, output format, merged output and ID stability.
Needs gffcompare and gffread on PATH (dev/environment.yml)."""

import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, '..', '..', 'bin')
SCRIPT = os.path.join(BIN, 'curate_models.py')
TARGET = 'Tgt_1'

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not (shutil.which('gffcompare') and shutil.which('gffread')), reason='needs gffcompare and gffread on PATH'
    ),
]


def model(tid, sp, chrom, strand, exons, de=0.1):
    tid = f'{tid}|{sp}'
    s, e = exons[0][0], exons[-1][1]
    rows = [
        f'{chrom}\tfomo\tgene\t{s}\t{e}\t1\t{strand}\t.\tID=gene-{tid};source={sp};gene_class=lncRNA',
        f'{chrom}\tfomo\ttranscript\t{s}\t{e}\t1\t{strand}\t.\tID={tid};Parent=gene-{tid};source={sp};gene_class=lncRNA;AS=100;de={de}',
    ]
    rows += [
        f'{chrom}\tfomo\texon\t{a}\t{b}\t1\t{strand}\t.\tID={tid}.exon{i};Parent={tid}'
        for i, (a, b) in enumerate(exons, 1)
    ]
    return rows


M = []
# 1  chain shared by A and B, different ends, equal exonic length -> rep by lowest de (b1)
M += model('a1', 'A', 'chr1', '+', [(100, 200), (300, 400), (500, 600)], de=0.2)
M += model('b1', 'B', 'chr1', '+', [(150, 200), (300, 400), (500, 650)], de=0.1)
# 2  same chain twice in one species -> unsupported
M += model('a2x', 'A', 'chr1', '+', [(1000, 1100), (1200, 1300)])
M += model('a2y', 'A', 'chr1', '+', [(1010, 1100), (1200, 1310)])
# 3  two supported chains sharing one SJ -> one gene, two transcripts: t1 = more SJs at equal species (a3)
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
# 8  R2 (3 sp) is t1, R1 (2 sp) t2, at equal SJ count; they share one SJ
M += model('a8r1', 'A', 'chr1', '+', [(5000, 5100), (5200, 5300), (5400, 5500), (5600, 5700)])
M += model('b8r1', 'B', 'chr1', '+', [(5000, 5100), (5200, 5300), (5400, 5500), (5600, 5700)])
M += model('a8r2', 'A', 'chr1', '+', [(5050, 5100), (5200, 5300), (5450, 5500), (5650, 5700)])
M += model('b8r2', 'B', 'chr1', '+', [(5050, 5100), (5200, 5300), (5450, 5500), (5650, 5700)])
M += model('c8r2', 'C', 'chr1', '+', [(5050, 5100), (5200, 5300), (5450, 5500), (5650, 5700)])
# 8b everything equal -> smallest ID; '-' strand: exon numbering 5'->3'
M += model('m_b', 'A', 'chr1', '-', [(6000, 6100), (6200, 6300)])
M += model('m_a', 'B', 'chr1', '-', [(6000, 6100), (6200, 6300)])
# 9  escaped comma in an ID survives as-is into supporting_transcripts
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
# 17 support first: F (1 SJ, 3 species) is t1, L (3 SJs, 2 species) t2; they share an SJ
M += model('a17l', 'A', 'chr1', '+', [(10000, 10100), (10200, 10300), (10400, 10500), (10600, 10700)])
M += model('b17l', 'B', 'chr1', '+', [(10000, 10100), (10200, 10300), (10400, 10500), (10600, 10700)])
M += model('c17f', 'C', 'chr1', '+', [(10050, 10100), (10200, 10350)])
M += model('d17f', 'D', 'chr1', '+', [(10050, 10100), (10200, 10350)])
M += model('e17f', 'E', 'chr1', '+', [(10050, 10100), (10200, 10350)])
# 18 ref_location classes on chr5 (reference transcripts in REF / REF_GTF)
LOC = {
    'a18so|A': 'sense_span_overlap',
    'a18ai|A': 'antisense_intronic',
    'a18dv|A': 'divergent',
    'a18dm|A': 'divergent',
    'a18far|A': 'intergenic',
    'a18cv|A': 'intergenic',
    'a6|A': 'antisense_exonic',
    'a7|A': 'intronic',
}
for sp in ('A', 'B'):
    lo = sp.lower()
    M += model(f'{lo}18so', sp, 'chr5', '+', [(1400, 1450), (1600, 1650)])  # straddles ref span 1500-1750 (+)
    M += model(f'{lo}18ai', sp, 'chr5', '-', [(3600, 3650), (3700, 3750)])  # inside ref (+) intron 3550-3800
    M += model(f'{lo}18dv', sp, 'chr5', '+', [(5000, 5100), (5200, 5300)])  # ref (-) ends 4800: 5' ends 200 bp apart
    M += model(f'{lo}18far', sp, 'chr5', '+', [(8000, 8100), (8200, 8300)])  # ref (-) ends 6990: 1010 bp > 1000
    M += model(f'{lo}18cv', sp, 'chr5', '+', [(13000, 13100), (13200, 13300)])  # ref (-) downstream: tail-to-tail
    M += model(f'{lo}18dm', sp, 'chr5', '-', [(20000, 20100), (20200, 20300)])  # ref (+) starts 20500: 200 bp
# 19 (plans/22) nearest reference gene tag on chr6; chr7 has no reference at all
TAG = {  # rep -> its gene's (ref_location, ref_gene_id, ref_gene_name, ref_gene_biotype, orientation, distance)
    'a61|A': ('intergenic', 'G1', 'GENEONE', 'protein_coding', 'sense', '700'),  # gene Name + gene_biotype
    'a62|A': ('intergenic', 'G2', 'REFSEQSYM', 'lnc_RNA', 'antisense', '400'),  # RefSeq gene=; biotype from type
    'a63|A': ('intergenic', 'G2', 'REFSEQSYM', 'lnc_RNA', 'sense', '1700'),  # tie 1700 vs G3 (+): sense wins
    'a64|A': ('intronic', 'G3', None, 'lncRNA', 'sense', '0'),  # exons parented by the gene; no name
    'a71|A': ('intergenic', None, None, None, None, None),  # seqid without reference: no ref_* keys
}
for sp in ('A', 'B'):
    lo = sp.lower()
    M += model(f'{lo}61', sp, 'chr6', '+', [(2000, 2100), (2200, 2300)])  # G1 (+) ends 1300 -> 700
    M += model(f'{lo}62', sp, 'chr6', '+', [(4300, 4400), (4500, 4600)])  # G2 (-) starts 5000 -> 400; G1 3000
    M += model(f'{lo}63', sp, 'chr6', '-', [(7000, 7100), (7200, 7300)])  # G2 (-) ends 5300 / G3 (+) starts 9000
    M += model(f'{lo}64', sp, 'chr6', '+', [(9120, 9130), (9160, 9180)])  # inside G3's intron 9101-9199
    M += model(f'{lo}71', sp, 'chr7', '+', [(100, 200), (300, 400)])
# 20 within a chain the longest model is the representative, even with a higher de
M += model('a20', 'A', 'chr1', '+', [(12000, 12100), (12200, 12500)], de=0.3)
M += model('b20', 'B', 'chr1', '+', [(12050, 12100), (12200, 12300)], de=0.01)
# 21 equal species and SJs -> the longer chain is t1 (Q); P and Q share their first SJ
M += model('a21p', 'A', 'chr1', '+', [(13000, 13100), (13200, 13300), (13400, 13500)])
M += model('b21p', 'B', 'chr1', '+', [(13000, 13100), (13200, 13300), (13400, 13500)])
M += model('a21q', 'A', 'chr1', '+', [(13000, 13100), (13200, 13350), (13450, 13800)])
M += model('b21q', 'B', 'chr1', '+', [(13000, 13100), (13200, 13350), (13450, 13800)])
# 22 gene location on the footprint: t1 (U, 3 sp) alone would be divergent; t2 (V) adds an exon
#    antisense to the chr4 (+) reference exon 1550-1560 -> the gene is antisense_exonic
for sp in ('A', 'B', 'C'):
    M += model(f'{sp.lower()}22u', sp, 'chr4', '-', [(1000, 1100), (1200, 1300)])
for sp in ('D', 'E'):
    M += model(f'{sp.lower()}22v', sp, 'chr4', '-', [(1000, 1100), (1200, 1300), (1500, 1600)])
LOC['a22u|A'] = 'antisense_exonic'
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
chr4\tRefSeq\tgene\t1550\t1560\t.\t+\t.\tID=gene-c4b;gene_biotype=protein_coding
chr4\tRefSeq\tmRNA\t1550\t1560\t.\t+\t.\tID=rna-c4b;Parent=gene-c4b
chr4\tRefSeq\texon\t1550\t1560\t.\t+\t.\tID=exon-c4b;Parent=rna-c4b
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
chr6\tRefSeq\tgene\t1000\t1300\t.\t+\t.\tID=G1;Name=GENEONE;gene_biotype=protein_coding
chr6\tRefSeq\tmRNA\t1000\t1300\t.\t+\t.\tID=t1;Parent=G1;Name=GENEONE-201
chr6\tRefSeq\texon\t1000\t1100\t.\t+\t.\tID=exon-t1-1;Parent=t1
chr6\tRefSeq\texon\t1200\t1300\t.\t+\t.\tID=exon-t1-2;Parent=t1
chr6\tRefSeq\tgene\t5000\t5300\t.\t-\t.\tID=G2
chr6\tRefSeq\tlnc_RNA\t5000\t5300\t.\t-\t.\tID=t2;Parent=G2;gene=REFSEQSYM
chr6\tRefSeq\texon\t5000\t5100\t.\t-\t.\tID=exon-t2-1;Parent=t2
chr6\tRefSeq\texon\t5200\t5300\t.\t-\t.\tID=exon-t2-2;Parent=t2
chr6\tRefSeq\tgene\t9000\t9300\t.\t+\t.\tID=G3;gene_biotype=lncRNA
chr6\tRefSeq\texon\t9000\t9100\t.\t+\t.\tID=exon-g3-1;Parent=G3
chr6\tRefSeq\texon\t9200\t9300\t.\t+\t.\tID=exon-g3-2;Parent=G3
"""

REF_GTF = (
    '\n'.join(
        '\t'.join(
            [c, 'ens', 'exon', s, e, '.', st, '.', f'gene_id "{t}"; transcript_id "{t}"; transcript_biotype "{bt}";']
        )
        for c, s, e, st, t, bt in [
            ('chr2', '700', '750', '+', 'trna', 'tRNA'),
            ('chr2', '1050', '1080', '+', 'm0', 'protein_coding'),
            ('chr2', '1900', '1950', '+', 'l7', 'lncRNA'),
            ('chr2', '2500', '2600', '+', 'l7', 'lncRNA'),
            ('chr2', '3050', '3060', '+', 'm1', 'protein_coding'),
            ('chr2', '3250', '3260', '+', 'as1', 'antisense'),
            ('chr2', '4050', '4060', '+', 'cds', 'protein_coding'),
            ('chr2', '5050', '5060', '+', 'misc', 'misc_RNA'),
            ('chr2', '6050', '6060', '+', 'mir', 'miRNA'),
            ('chr3', '5000', '6000', '+', 'c3', 'protein_coding'),
            ('chr4', '1550', '1560', '+', 'c4b', 'protein_coding'),
            ('chr5', '1500', '1550', '+', 'so', 'lncRNA'),
            ('chr5', '1700', '1750', '+', 'so', 'lncRNA'),
            ('chr5', '3500', '3550', '+', 'ai', 'protein_coding'),
            ('chr5', '3800', '3850', '+', 'ai', 'protein_coding'),
            ('chr5', '4500', '4600', '-', 'dv', 'protein_coding'),
            ('chr5', '4700', '4800', '-', 'dv', 'protein_coding'),
            ('chr5', '6500', '6600', '-', 'far', 'lncRNA'),
            ('chr5', '6800', '6990', '-', 'far', 'lncRNA'),
            ('chr5', '13500', '13600', '-', 'cv', 'protein_coding'),
            ('chr5', '13700', '13800', '-', 'cv', 'protein_coding'),
            ('chr5', '20500', '20600', '+', 'dm', 'protein_coding'),
            ('chr5', '20700', '20800', '+', 'dm', 'protein_coding'),
        ]
    )
    + '\n'
    + ''.join(
        '\t'.join([c, 'ens', 'exon', s, e, '.', st, '.', att]) + '\n'
        for c, s, e, st, att in [
            (
                'chr6',
                '1000',
                '1100',
                '+',
                'gene_id "G1"; transcript_id "t1"; gene_name "GENEONE"; gene_biotype "protein_coding";',
            ),
            (
                'chr6',
                '1200',
                '1300',
                '+',
                'gene_id "G1"; transcript_id "t1"; gene_name "GENEONE"; gene_biotype "protein_coding";',
            ),
            ('chr6', '5000', '5100', '-', 'gene_id "G2"; transcript_id "t2"; gene_name "REFSEQSYM";'),
            ('chr6', '5200', '5300', '-', 'gene_id "G2"; transcript_id "t2"; gene_name "REFSEQSYM";'),
            ('chr6', '9000', '9100', '+', 'gene_id "G3"; transcript_id "t3"; gene_type "lncRNA";'),
            ('chr6', '9200', '9300', '+', 'gene_id "G3"; transcript_id "t3"; gene_type "lncRNA";'),
        ]
    )
)


def run(tmp, name, gff_rows, extra, ref=None):
    gff = os.path.join(tmp, f'{TARGET}.allModels.lnc_RNA.raw.gff3')
    with open(gff, 'w') as o:
        o.write('##gff-version 3\n' + ''.join(r + '\n' for r in gff_rows))
    prefix = os.path.join(tmp, name, f'{TARGET}.curated.lnc_RNA')
    cmd = [sys.executable, SCRIPT, '--gff', gff, '--prefix', prefix] + extra + (['--ref', ref] if ref else [])
    return subprocess.run(cmd, capture_output=True, text=True), prefix


def read_gff(path):
    """-> genes {ID: (fields, attrs)}, transcripts {representative model: (fields, attrs)}, exon rows
    [(fields, attrs)], in file order."""
    genes, tx, exons = {}, {}, []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line[0] == '#':
                continue
            f = line.rstrip('\n').split('\t')
            a = dict(kv.split('=', 1) for kv in f[8].split(';'))
            if f[2] == 'gene':
                genes[a['ID']] = (f, a)
            elif f[2] == 'exon':
                exons.append((f, a))
            else:
                tx[a['representative_transcript']] = (f, a)
    return genes, tx, exons


def gene_of(genes, tx, rid):
    """-> attrs of the gene whose transcript is represented by model rid ({} if none)."""
    t = tx.get(rid)
    return genes[t[1]['Parent']][1] if t else {}


def read_text(path):
    with gzip.open(path, 'rt') as fh:
        return fh.read()


@pytest.fixture(scope='module')
def tmp(tmp_path_factory):
    return str(tmp_path_factory.mktemp('curate'))


@pytest.fixture(scope='module')
def ref(tmp):
    path = os.path.join(tmp, 'ref.gff3')
    with open(path, 'w') as o:
        o.write(REF)
    return path


@pytest.fixture(scope='module')
def ref_gtf(tmp):
    path = os.path.join(tmp, 'ref.gtf')
    with open(path, 'w') as o:
        o.write(REF_GTF)
    return path


@pytest.fixture(scope='module')
def main(tmp, ref):
    """The main run: every fixture model, GFF3 reference."""
    r, pre = run(tmp, 'main', M, [], ref)
    assert r.returncode == 0, r.stderr
    genes, tx, exons = read_gff(pre + '.gff3.gz')
    with open(pre + '.curation.json') as fh:
        report = json.load(fh)
    return SimpleNamespace(
        pre=pre,
        genes=genes,
        tx=tx,
        exons=exons,
        rep=set(tx),
        o=report['observed'],
        d=report['discarded_by_reference'],
        n=report['novel'],
        G=lambda rid: gene_of(genes, tx, rid),
        T=lambda rid: tx.get(rid, (None, {}))[1],
    )


@pytest.fixture(scope='module')
def cm():
    sys.path.insert(0, BIN)
    import curate_models

    return curate_models


EXPECTED_REPRESENTATIVES = (
    {
        'b1|B',
        'a3|A',
        'c3|C',
        'a4p|A',
        'c4q|C',
        'a8r1|A',
        'a8r2|A',
        'm_a|B',
        'x%2Cy|B',
        'a6|A',
        'a7|A',
        'a15|A',
        'a17l|A',
        'c17f|C',
        'a20|A',
        'a21p|A',
        'a21q|A',
        'd22v|D',
    }
    | set(LOC)
    | set(TAG)
)


def is_tn(attrs, n):
    return attrs.get('ID', '').startswith(f'fomo_lncRNA_t{n}_')


def test_chain_support(main):
    rep, tx_of, o = main.rep, main.T, main.o
    assert 'b1|B' in rep and 'a1|A' not in rep, '1  shared chain kept; equal length -> lowest de (b1)'
    assert not any(x.startswith(('a2x', 'a2y')) for x in rep), '2  same-species duplicate chain is not supported'
    assert tx_of('a3|A').get('Parent') and tx_of('a3|A').get('Parent') == tx_of('c3|C').get('Parent'), (
        '3  chains sharing one SJ -> one gene with two transcripts (a3, c3)'
    )
    assert len([x for x in rep if x[:3] in ('a4p', 'b4p', 'c4q', 'd4q')]) == 2 and tx_of('a4p|A').get(
        'Parent'
    ) != tx_of('c4q|C').get('Parent'), '4  unsupported bridge does not merge two genes'
    assert tx_of('x%2Cy|B').get('supporting_transcripts') == 'x%2Cy|B,xy|C', (
        '9  escaped comma in ID kept verbatim in supporting_transcripts'
    )
    assert o['single_exon'] == 2, '10 single-exon models counted and skipped'
    assert o['self_excluded'] == 1 and not any(x.startswith(('t14', 'a14')) for x in rep), (
        '14 self projection gives no support'
    )


def test_transcript_ranking(main):
    tx_of, gene_attrs = main.T, main.G
    assert is_tn(tx_of('a3|A'), 1) and is_tn(tx_of('c3|C'), 2), '3  equal species -> more SJs is t1'
    assert gene_attrs('a3|A').get('n_chains') == '2' and gene_attrs('a3|A').get('n_species') == '4', (
        '3  gene reports n_chains=2, n_species=4'
    )
    assert (
        tx_of('a3|A').get('supporting_species') == 'A,B'
        and tx_of('a3|A').get('supporting_transcripts') == 'a3|A,b3|B'
        and tx_of('c3|C').get('supporting_transcripts') == 'c3|C,d3|D'
    ), '3  each transcript lists its own chain support'
    assert (
        is_tn(tx_of('a8r2|A'), 1)
        and is_tn(tx_of('a8r1|A'), 2)
        and tx_of('a8r1|A').get('Parent') == tx_of('a8r2|A').get('Parent')
    ), '8  equal SJ count -> more supporting species is t1'
    assert 'm_a|B' in main.rep and 'm_b|A' not in main.rep, '8b full tie -> smallest ID'
    assert is_tn(tx_of('c17f|C'), 1) and is_tn(tx_of('a17l|A'), 2), (
        '17 support first: F (1 SJ, 3 sp) is t1, L (3 SJs, 2 sp) is t2'
    )
    assert (
        tx_of('a17l|A').get('n_sj') == '3'
        and tx_of('c17f|C').get('n_sj') == '1'
        and gene_attrs('a17l|A').get('n_chains') == '2'
        and gene_attrs('a17l|A').get('n_species') == '5'
        and tx_of('c17f|C').get('supporting_species') == 'C,D,E'
    ), '17 transcripts report their own SJs / species; gene n_chains=2, n_species=5'
    assert 'a20|A' in main.rep and 'b20|B' not in main.rep, '20 longest model is the representative despite a higher de'
    assert is_tn(tx_of('a21q|A'), 1) and is_tn(tx_of('a21p|A'), 2), '21 equal species and SJs -> the longer chain is t1'
    assert is_tn(tx_of('a22u|A'), 1) and tx_of('d22v|D').get('Parent') == tx_of('a22u|A').get('Parent'), (
        '22 U (3 sp) is t1, V in the same gene'
    )
    g22 = gene_attrs('a22u|A')
    assert g22.get('ref_location') == 'antisense_exonic' and g22.get('ref_gene_id') == 'gene-c4b', (
        "22 gene location is computed on the footprint (t2's antisense exon), not on t1"
    )
    gf = main.genes.get(tx_of('a22u|A').get('Parent'), ([None] * 9,))[0]
    assert gf[3] == '1000' and gf[4] == '1600', '22 gene span covers all its transcripts'


def test_reference_filter(main):
    rep, d, n = main.rep, main.d, main.n
    assert not any(x.startswith(('a5', 'b5')) for x in rep), '5  ref exon removes b5 -> chain loses support'
    assert d.get('chains_lost_by_type:tRNA') == 1, '5  lost chain reported under tRNA'
    assert 'a6|A' in rep or 'b6|B' in rep, '6  opposite-strand ref exon does not remove (same-strand rule)'
    assert n.get('genes_location:antisense_exonic') == 2, '6  counted as antisense_exonic (with 22)'
    assert 'a7|A' in rep or 'b7|B' in rep, '7  inside a ref intron is kept'
    assert n.get('genes_location:intronic') == 2, '7  counted as intronic (with 19 a64)'
    assert d.get('chains_lost_by_type:mRNA') == 2 and d.get('chains_lost_by_type:antisense_RNA') == 1, (
        '16 mRNA + antisense_RNA chain counted once per type (mRNA also from the CDS-only model)'
    )
    assert d.get('chains_lost_by_primary_type:mRNA') == 2 and 'chains_lost_by_primary_type:antisense_RNA' not in d, (
        '16 primary type: mRNA wins over antisense_RNA'
    )
    assert d.get('chains_lost_by_type:transcript:misc_RNA') == 1, (
        '16 RefSeq misc_RNA transcript labelled transcript:misc_RNA'
    )
    assert d.get('chains_lost_by_type:miRNA') == 1, '16 childless miRNA counted by its span'
    tot_primary = sum(v for k, v in d.items() if k.startswith('chains_lost_by_primary_type:'))
    assert tot_primary == d['supported_chains_lost'], '16 primary-type counts add up to supported_chains_lost'


def test_ref_location_and_tags(main):
    gene_attrs, genes, tx, n = main.G, main.genes, main.tx, main.n
    for rid, want in LOC.items():
        assert gene_attrs(rid).get('ref_location') == want, (
            f"18 {rid} ref_location={want} (got {gene_attrs(rid).get('ref_location')})"
        )
    assert not any(k.startswith('ref_') for _, a in tx.values() for k in a), '18 ref_* attributes only on gene rows'
    far = gene_attrs('a18far|A')
    assert (
        far.get('ref_distance') == '1010'
        and far.get('ref_gene_id') == 'rna-far'
        and far.get('ref_gene_biotype') == 'lnc_RNA'
        and far.get('ref_gene_orientation') == 'antisense'
    ), '18 intergenic just beyond --divergent-dist: nearest rna-far (lnc_RNA, antisense), ref_distance=1010'
    for rid, want in TAG.items():
        a = gene_attrs(rid)
        got = tuple(
            a.get(k)
            for k in (
                'ref_location',
                'ref_gene_id',
                'ref_gene_name',
                'ref_gene_biotype',
                'ref_gene_orientation',
                'ref_distance',
            )
        )
        assert got == want, f'19 {rid} tag {got}'
    assert not any('ref_nearest' in a for _, a in genes.values()), '19 ref_nearest is not written'
    assert not any(v == '.' for _, a in list(genes.values()) + list(tx.values()) for v in a.values()), (
        "19 no attribute is written as '.'"
    )
    inter = [
        int(a['ref_distance'])
        for _, a in genes.values()
        if a.get('ref_location') == 'intergenic' and 'ref_distance' in a
    ]
    assert (
        n.get('intergenic_distance_n') == len(inter)
        and n.get('intergenic_distance_min') == min(inter)
        and n.get('intergenic_distance_max') == max(inter)
    ), '19 intergenic distance stats'
    dv = gene_attrs('a18dv|A')
    assert dv.get('ref_distance') == '200', '18 divergent: ref_distance=200'
    assert dv.get('ref_gene_id') == 'rna-dv' and dv.get('ref_gene_orientation') == 'antisense', (
        '18 divergent: tagged with its head-to-head partner rna-dv (antisense)'
    )
    assert all(
        a.get('ref_distance') == '0'
        for _, a in genes.values()
        if a.get('ref_location') in ('antisense_exonic', 'intronic', 'sense_span_overlap', 'antisense_intronic')
    ), '18 overlapping classes have ref_distance=0'
    loc_tot = sum(v for k, v in n.items() if k.startswith('genes_location:'))
    nloc = {k.split(':', 1)[1]: v for k, v in n.items() if k.startswith('genes_location:')}
    assert loc_tot == n['genes'] and len(nloc) == 7 and nloc['no_reference'] == 0 and nloc['divergent'] == 2, (
        f'18 all 7 genes_location classes written, sum to genes ({nloc})'
    )


def test_expected_representatives_and_report(main):
    assert main.rep == EXPECTED_REPRESENTATIVES, f'{sorted(main.rep ^ EXPECTED_REPRESENTATIVES)}'
    n = main.n
    assert (
        n['genes'] == len(main.genes)
        and n['transcripts'] == len(main.tx) == len(EXPECTED_REPRESENTATIVES)
        and n['genes'] == len(EXPECTED_REPRESENTATIVES) - 5
    )


def test_output_format(main):
    genes, tx, exons, tx_of = main.genes, main.tx, main.exons, main.T
    assert all(a.get('biotype') == 'lncRNA' for _, a in list(genes.values()) + list(tx.values())), (
        'biotype=lncRNA on genes and transcripts'
    )
    assert all(f[2] == 'lnc_RNA' and f[5] == '.' for f, _ in tx.values()) and all(
        f[5] == '.' for f, _ in genes.values()
    ), 'transcripts typed lnc_RNA; column 6 is . everywhere'
    assert all(
        set(a)
        == {
            'ID',
            'Parent',
            'biotype',
            'representative_transcript',
            'supporting_transcripts',
            'supporting_species',
            'n_sj',
        }
        for _, a in tx.values()
    ), 'transcript attributes exactly the set'
    assert all(
        a['ID'] == f"fomo_lncRNA_g_{a['ID'][-12:]}" and len(a['ID']) == len('fomo_lncRNA_g_') + 12
        for _, a in genes.values()
    ), 'gene IDs fomo_lncRNA_g_<H>'
    assert all(
        a['ID'] == f"fomo_lncRNA_t{a['ID'].split('_t')[1].split('_')[0]}_{a['Parent'][-12:]}" for _, a in tx.values()
    ), 'transcript IDs fomo_lncRNA_t<n>_<H> share their gene hash'
    for f, a in tx.values():
        rows = [(ef, ea) for ef, ea in exons if ea['Parent'] == a['ID']]
        nums = [int(re.search(r'_t\d+e(\d+)_', ea['ID']).group(1)) for _, ea in rows]
        want = list(range(1, len(rows) + 1)) if f[6] == '+' else list(range(len(rows), 0, -1))
        assert nums == want, f"exon numbering of {a['ID']}: {nums}"
        assert [int(ef[3]) for ef, _ in rows] == sorted(int(ef[3]) for ef, _ in rows), f"exon order of {a['ID']}"
    me = [ea['ID'] for _, ea in exons if ea['Parent'] == tx_of('m_a|B').get('ID')]
    assert me and me[0].startswith('fomo_lncRNA_t1e2_') and me[1].startswith('fomo_lncRNA_t1e1_'), (
        "8b '-' strand: rows t1e2, t1e1"
    )
    ids = list(genes) + [a['ID'] for _, a in tx.values()] + [a['ID'] for _, a in exons]
    assert len(ids) == len(set(ids)), 'every ID unique'


def test_gffread_reads_curated(main, tmp):
    plain = os.path.join(tmp, 'curated_plain.gff3')
    with open(plain, 'w') as o:
        o.write(read_text(main.pre + '.gff3.gz'))
    g = subprocess.run(['gffread', '-E', '-T', '-o', '-', plain], capture_output=True, text=True)
    gtx = {ln.split('transcript_id "')[1].split('"')[0] for ln in g.stdout.splitlines() if 'transcript_id' in ln}
    assert g.returncode == 0 and gtx == {a['ID'] for _, a in main.tx.values()}, 'one transcript per chain'


def test_merged_file(main):
    merged = main.pre + '.merged.gff3.gz'
    assert os.path.exists(merged), '15 merged file written'
    lines = read_text(merged).splitlines()
    feat = [ln.split('\t') for ln in lines if ln and ln[0] != '#']
    ids = [kv[3:] for f in feat for kv in f[8].split(';') if kv.startswith('ID=')]
    assert len(ids) == len(set(ids)), '15 merged IDs unique'
    order = [f[0] for f in feat]
    blocks = [c for i, c in enumerate(order) if i == 0 or c != order[i - 1]]
    assert len(blocks) == len(set(blocks)), f'15 merged seqids contiguous ({blocks})'
    assert blocks == ['chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr1', 'chr7'], (
        f'15 curated chr2/3/4/5/6 genes inside their seqid block, chr1/chr7 appended ({blocks})'
    )
    ref_lines = [ln for ln in REF.splitlines() if ln]
    assert [ln for ln in lines if '\tfomo\t' not in ln and not ln.startswith('# curate_models')] == ref_lines, (
        '15 reference lines verbatim and in order'
    )
    r2 = subprocess.run(
        ['gffread', '-E', '-o', os.devnull, '/dev/stdin'], input='\n'.join(lines) + '\n', capture_output=True, text=True
    )
    assert r2.returncode == 0, '15 gffread accepts merged'


def test_merged_without_gff_version_header(tmp):
    ref_nh = os.path.join(tmp, 'ref_noheader.gff3')
    with open(ref_nh, 'w') as o:
        o.write(REF.split('\n', 1)[1])
    r, pre = run(tmp, 'noheader', M, [], ref_nh)
    lines = read_text(pre + '.merged.gff3.gz').splitlines()
    ref_lines = [ln for ln in REF.splitlines() if ln]
    assert (
        r.returncode == 0
        and lines[0] == '##gff-version 3'
        and lines[1].startswith('# curate_models')
        and [ln for ln in lines[2:] if '\tfomo\t' not in ln] == ref_lines[1:]
    ), '15 reference without ##gff-version: header added, reference lines verbatim'


def test_no_reference(tmp):
    r, pre = run(tmp, 'noref', M, [])
    assert r.returncode == 0 and not os.path.exists(pre + '.merged.gff3.gz'), '12 no reference -> no merged file'
    genes, tx, _ = read_gff(pre + '.gff3.gz')
    assert set(tx) == EXPECTED_REPRESENTATIVES | {'b5|B', 'a16|A', 'a16c|A', 'a16m|A', 'a16r|A'}, (
        '12 without reference every supported chain is kept'
    )
    assert all(
        a.get('ref_location') == 'no_reference' and not any(k.startswith('ref_') and k != 'ref_location' for k in a)
        for _, a in genes.values()
    ), '12 without reference: ref_location=no_reference, no other ref_* key'
    with open(pre + '.curation.json') as fh:
        novel = json.load(fh)['novel']
    assert novel.get('genes_location:no_reference') == novel['genes'] and novel.get('genes_location:intergenic') == 0, (
        '12 every location class written; all genes no_reference'
    )


def test_keep_self(tmp):
    _, pre = run(tmp, 'keepself', M, ['--keep-self'])
    _, tx, _ = read_gff(pre + '.gff3.gz')
    assert any(x in ('t14|Tgt_1', 'a14|A') for x in tx), '14 --keep-self lets the target support'


def test_gtf_reference(tmp, ref_gtf, main):
    r, pre = run(tmp, 'gtf', M, [], ref_gtf)
    assert r.returncode == 0, '13 GTF reference run exits 0'
    genes, tx, _ = read_gff(pre + '.gff3.gz')
    assert set(tx) == EXPECTED_REPRESENTATIVES, '13 GTF reference gives the same curated set'

    def g(rid):
        return gene_of(genes, tx, rid)

    assert all(g(k).get('ref_location') == v for k, v in LOC.items()), (
        '13 GTF reference gives the same ref_location classes'
    )
    assert (
        g('a61|A').get('ref_gene_name') == 'GENEONE'
        and g('a61|A').get('ref_gene_biotype') == 'protein_coding'
        and g('a62|A').get('ref_gene_name') == 'REFSEQSYM'
        and g('a64|A').get('ref_gene_biotype') == 'lncRNA'
        and g('a64|A').get('ref_gene_id') == 'G3'
        and 'ref_gene_id' not in g('a71|A')
    ), '19 GTF reference: gene_name / gene_biotype / gene_type and gene_id tags'
    assert set(genes) == set(main.genes), '13 IDs do not depend on the reference format'


@pytest.mark.parametrize(
    'loc, bad',
    [
        ('intronic', dict(gene='g', name='.', biotype='.', orientation='sense', distance=5)),
        ('antisense_exonic', dict(gene='g', name='.', biotype='.', orientation='sense', distance=0)),
        ('divergent', dict(gene='g', name='.', biotype='.', orientation='antisense', distance=5000)),
    ],
)
def test_location_guard_rejects_contradictions(cm, loc, bad):
    with pytest.raises(cm.fomo_io.FomoError, match='inconsistent reference tag'):
        cm.check_location(loc, bad, 1000, 'x')


def test_empty_input(tmp, ref):
    r, pre = run(tmp, 'empty', [], [], ref)
    assert r.returncode == 0, '11 empty input exits 0'
    genes, _, _ = read_gff(pre + '.gff3.gz')
    merged = [ln for ln in read_text(pre + '.merged.gff3.gz').splitlines() if not ln.startswith('# curate_models')]
    assert not genes and merged == [ln for ln in REF.splitlines() if ln], (
        '11 empty input -> no genes, merged == reference'
    )


def test_track_tags(tmp):
    _, pre6 = run(tmp, 'decoy', M[:10], ['--track', 'lnc_RNA.decoy'])
    g6, tx6, _ = read_gff(pre6 + '.gff3.gz')
    assert (
        g6
        and all(a['biotype'] == 'decoy_lncRNA' and a['ID'].startswith('fomo_decoy_lncRNA_g_') for _, a in g6.values())
        and all(f[2] == 'lnc_RNA' for f, _ in tx6.values())
    ), 'decoy: biotype decoy_lncRNA, IDs fomo_decoy_lncRNA_*, type lnc_RNA'
    _, pre7 = run(tmp, 'mrna', M[:10], ['--track', 'mRNA'])
    g7, tx7, _ = read_gff(pre7 + '.gff3.gz')
    assert (
        g7
        and all(a['biotype'] == 'mRNA' and a['ID'].startswith('fomo_mRNA_g_') for _, a in g7.values())
        and all(f[2] == 'mRNA' for f, _ in tx7.values())
    ), 'mRNA: biotype mRNA, IDs fomo_mRNA_*, type mRNA'
    assert set(g6).isdisjoint(g7) and {k[-12:] for k in g6}.isdisjoint(k[-12:] for k in g7), (
        'the same chains in two tracks get different hashes'
    )


def test_id_clash_with_reference(tmp, main):
    ref_c = os.path.join(tmp, 'ref_clash.gff3')
    with open(ref_c, 'w') as o:
        o.write(REF + f'chr4\tRefSeq\tgene\t300\t400\t.\t+\t.\tID={sorted(main.genes)[0]}\n')
    r, _ = run(tmp, 'clash', M, [], ref_c)
    assert r.returncode == 1 and 'already used in the reference' in r.stderr, (
        'curated ID clashing with the reference exits 1'
    )
    r, pre = run(tmp, 'clash2', M, ['--id-prefix', 'X'], ref_c)
    genes, _, _ = read_gff(pre + '.gff3.gz')
    assert r.returncode == 0 and all(k.startswith('X_g_') for k in genes), '--id-prefix X avoids the clash'


def test_duplicate_transcript_id_exits_1(tmp):
    r, _ = run(tmp, 'dup', M[:6] + [M[1]], [])
    assert r.returncode == 1 and 'duplicate transcript ID' in r.stderr


def test_rerun_is_byte_identical(tmp, ref, main):
    _, pre = run(tmp, 'main2', M, [], ref)
    for suffix in ('.gff3.gz', '.merged.gff3.gz'):
        assert read_text(main.pre + suffix) == read_text(pre + suffix), suffix


def test_ids_stable_across_input_order(tmp, ref, main):
    models = []
    for ln in M:
        if '\tgene\t' in ln:
            models.append([])
        models[-1].append(ln)
    shuffled = [ln for m in reversed(models) for ln in m]
    r, pre = run(tmp, 'shuffled', shuffled, [], ref)
    _, tx, _ = read_gff(pre + '.gff3.gz')
    assert r.returncode == 0 and {k: a['ID'] for k, (_, a) in tx.items()} == {
        k: a['ID'] for k, (_, a) in main.tx.items()
    }


def test_ids_stable_with_more_sources(tmp, ref, main):
    extra = (
        model('f1', 'F', 'chr1', '+', [(120, 200), (300, 400), (500, 610)])
        + model('a90', 'A', 'chr8', '+', [(100, 200), (300, 400)])
        + model('b90', 'B', 'chr8', '+', [(100, 200), (300, 400)])
    )
    r, pre = run(tmp, 'moresources', M + extra, [], ref)
    genes, tx, _ = read_gff(pre + '.gff3.gz')
    assert (
        r.returncode == 0
        and set(main.genes) < set(genes)
        and main.T('b1|B')['ID'] == tx['b1|B'][1]['ID']
        and 'F' in tx['b1|B'][1]['supporting_species']
    )


def test_gene_hash_collision_exits_1(tmp, cm, monkeypatch):
    gff = os.path.join(tmp, f'{TARGET}.allModels.lnc_RNA.raw.gff3')
    with open(gff, 'w') as o:
        o.write('##gff-version 3\n' + ''.join(r + '\n' for r in M))
    buf = io.StringIO()
    monkeypatch.setattr(cm, 'gene_hash', lambda *a, **k: 'deadbeef0000')
    monkeypatch.setattr(
        sys,
        'argv',
        ['curate_models.py', '--gff', gff, '--prefix', os.path.join(tmp, 'collision', 'c'), '--skip-gffread'],
    )
    monkeypatch.setattr(sys, 'stderr', buf)
    with pytest.raises(SystemExit) as ex:
        cm.fomo_io.run(cm.main)
    assert ex.value.code == 1 and 'gene hash collision deadbeef0000' in buf.getvalue()
