"""bin/td2_coding_filter.py: drop lncRNA records with a complete, PSAURON-confident ORF."""

import pytest

import td2_coding_filter as td2
from conftest import run_script, write

PEP = """>t1|lnc_RNA|S.p1 GENE.t1~~t1.p1  ORF type:complete (+),psauron_score=0.91 len:120
MKV*
>t1|lnc_RNA|S.p2 GENE.t1~~t1.p2  ORF type:complete (+),psauron_score=0.20 len:100
MKV*
>t2|lnc_RNA|S.p1 GENE.t2~~t2.p1  ORF type:5prime_partial (+),psauron_score=0.99 len:100
MKV*
>t3|lnc_RNA|S.p1 GENE.t3~~t3.p1  ORF type:complete (+),psauron_score=0.50 len:100
MKV*
>gone|lnc_RNA|S.p1 GENE.x~~x.p1  ORF type:complete (+),psauron_score=0.99 len:100
MKV*
"""
FASTA = '>t1|lnc_RNA|S desc\nACGT\nACGT\n>t2|lnc_RNA|S\nAAAA\n>t3|lnc_RNA|S\nCCCC\n>t4|lnc_RNA|S\nGGGG\n'


def test_parse_coding_thresholds(tmp_path):
    pep = write(tmp_path / 'p.pep', PEP)
    assert td2.parse_coding(str(pep), 0.5) == {'t1|lnc_RNA|S', 't3|lnc_RNA|S', 'gone|lnc_RNA|S'}
    assert td2.parse_coding(str(pep), 0.95) == {'gone|lnc_RNA|S'}
    assert td2.parse_coding(str(pep), 1.01) == set(), 'a cut-off above 1 disables the filter'


@pytest.mark.parametrize(
    'stage, coding_ids', [('input', ['t1', 't3']), ('projected', ['t1|lnc_RNA|S', 't3|lnc_RNA|S'])]
)
def test_cli(tmp_path, stage, coding_ids):
    fa, pep = write(tmp_path / 'in.fa', FASTA), write(tmp_path / 'p.pep', PEP)
    kept, ids, mqc = tmp_path / 'kept.fa', tmp_path / 'ids.txt', tmp_path / 'x_td2_coding_mqc.tsv'
    r = run_script(
        'td2_coding_filter.py',
        '--in-fasta',
        fa,
        '--pep',
        pep,
        '--kept-fasta',
        kept,
        '--coding-ids',
        ids,
        '--mqc',
        mqc,
        '--sample',
        'S',
        '--stage',
        stage,
    )
    assert r.returncode == 0, r.stderr
    assert kept.read_text() == '>t2|lnc_RNA|S\nAAAA\n>t4|lnc_RNA|S\nGGGG\n', 'kept records verbatim, in order'
    assert ids.read_text().split() == coding_ids, 'input stage strips to the bare tid, projected keeps the id'
    assert mqc.read_text().splitlines()[1] == f'S\t{stage}\tlncRNA\t4\t2\t2\t50.00', 'hits outside the FASTA ignored'


def test_empty_pep_keeps_everything(tmp_path):
    fa, pep = write(tmp_path / 'in.fa', FASTA), write(tmp_path / 'p.pep', '')
    kept = tmp_path / 'kept.fa'
    run_script(
        'td2_coding_filter.py',
        '--in-fasta',
        fa,
        '--pep',
        pep,
        '--kept-fasta',
        kept,
        '--coding-ids',
        tmp_path / 'i',
        '--mqc',
        tmp_path / 'm',
        '--sample',
        'S',
        '--stage',
        'input',
        check=True,
    )
    assert kept.read_text() == FASTA
