"""bin/bam_to_gff.sh: spliced BAM alignments -> projected GFF3 models with globally unique IDs."""

import subprocess

import pytest

from conftest import needs, run_script, write

pytestmark = [pytest.mark.integration, needs('samtools')]

HEADER = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n'


def sam_record(qname, pos, cigar, flag=0, tags='NM:i:2\tAS:i:150\tde:f:0.01'):
    seq_len = sum(int(n) for n, op in _cigar(cigar) if op in 'MIS=X')
    return '\t'.join([qname, str(flag), 'chr1', str(pos), '60', cigar, '*', '0', '0', 'A' * seq_len, '*', tags])


def _cigar(cigar):
    import re

    return re.findall(r'(\d+)([MIDNSHP=X])', cigar)


def bam(tmp_path, records):
    sam = write(tmp_path / 'a.sam', HEADER + '\n'.join(records) + '\n')
    out = tmp_path / 'a.bam'
    subprocess.run(['samtools', 'view', '-b', '-o', str(out), str(sam)], check=True)
    return out


def rows(text, ftype=None):
    out = [ln.split('\t') for ln in text.splitlines() if ln and not ln.startswith('#')]
    return [r for r in out if ftype is None or r[2] == ftype]


def test_cigar_to_exons_and_ids(tmp_path):
    b = bam(tmp_path, [sam_record('t1|lnc_RNA|SpeciesA', 101, '5S50M2D10M1000N40M3I20M')])
    r = run_script('bam_to_gff.sh', b)
    assert r.returncode == 0, r.stderr
    assert r.stdout.startswith('##gff-version 3\n')
    exons = rows(r.stdout, 'exon')
    assert [(int(e[3]), int(e[4])) for e in exons] == [(101, 162), (1163, 1222)], (
        'D extends the block, N splits it, I adds nothing'
    )
    t = rows(r.stdout, 'transcript')[0]
    attrs = dict(kv.split('=', 1) for kv in t[8].split(';'))
    assert attrs['ID'] == 't1|SpeciesA' and attrs['Parent'] == 'gene-t1|SpeciesA'
    assert attrs['source'] == 'SpeciesA' and attrs['gene_class'] == 'lnc_RNA'
    assert attrs['NM'] == '2' and attrs['AS'] == '150' and attrs['de'] == '0.01'
    assert attrs['query_length'] == '128' and attrs['aligned_bases'] == '120'


def test_single_exon_dropped_unless_requested(tmp_path):
    b = bam(tmp_path, [sam_record('t1|lnc_RNA|A', 101, '100M'), sam_record('t2|lnc_RNA|A', 501, '50M100N50M')])
    default = run_script('bam_to_gff.sh', b)
    assert [r[8].split(';')[0] for r in rows(default.stdout, 'transcript')] == ['ID=t2|A']
    keep = run_script('bam_to_gff.sh', '--include-single-exon', b)
    assert len(rows(keep.stdout, 'transcript')) == 2


def test_secondary_supplementary_unmapped_skipped(tmp_path):
    b = bam(
        tmp_path,
        [
            sam_record('t1|lnc_RNA|A', 101, '50M100N50M'),
            sam_record('t1|lnc_RNA|A', 501, '50M100N50M', flag=256),
            sam_record('t1|lnc_RNA|A', 901, '50M100N50M', flag=2048),
        ],
    )
    r = run_script('bam_to_gff.sh', b)
    assert len(rows(r.stdout, 'transcript')) == 1


def test_same_tid_from_two_sources_stays_unique(tmp_path):
    b = bam(
        tmp_path,
        [
            sam_record('lnc_RNA1412|lnc_RNA|A', 101, '50M100N50M'),
            sam_record('lnc_RNA1412|lnc_RNA|B', 101, '50M100N50M'),
        ],
    )
    ids = [r[8].split(';')[0] for r in rows(run_script('bam_to_gff.sh', b).stdout, 'transcript')]
    assert ids == ['ID=lnc_RNA1412|A', 'ID=lnc_RNA1412|B']


def test_minus_strand(tmp_path):
    b = bam(tmp_path, [sam_record('t1|lnc_RNA|A', 101, '50M100N50M', flag=16)])
    assert {r[6] for r in rows(run_script('bam_to_gff.sh', b).stdout)} == {'-'}


def test_tid_with_gff_separator_fails(tmp_path):
    b = bam(tmp_path, [sam_record('t;1|lnc_RNA|A', 101, '50M100N50M')])
    r = run_script('bam_to_gff.sh', b)
    assert r.returncode != 0 and 'corrupt GFF attributes' in r.stderr
