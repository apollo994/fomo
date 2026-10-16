"""bin/td2_coding_filter.py: drop lncRNA records with a complete, PSAURON-confident ORF."""

import json

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


def run(tmp_path, fasta, pep, stage='input'):
    prefix = tmp_path / 'x'
    identity = '{"species": "S", "side": "source", "track": "lnc_RNA"}'
    args = ['--fasta', fasta, '--pep', pep, '--prefix', prefix, '--stage', stage, '--identity', identity]
    r = run_script('td2_coding_filter.py', *args)
    assert r.returncode == 0, r.stderr
    return r, prefix


def test_coding_names_thresholds(tmp_path):
    pep = write(tmp_path / 'p.pep', PEP)
    assert td2.coding_names(str(pep), 0.5) == {'t1|lnc_RNA|S', 't3|lnc_RNA|S', 'gone|lnc_RNA|S'}
    assert td2.coding_names(str(pep), 0.95) == {'gone|lnc_RNA|S'}
    assert td2.coding_names(str(pep), 1.01) == set(), 'a cut-off above 1 disables the filter'


@pytest.mark.parametrize(
    'stage, coding_ids', [('input', ['t1', 't3']), ('projected', ['t1|lnc_RNA|S', 't3|lnc_RNA|S'])]
)
def test_cli(tmp_path, stage, coding_ids):
    fa, pep = write(tmp_path / 'in.fa', FASTA), write(tmp_path / 'p.pep', PEP)
    _, prefix = run(tmp_path, fa, pep, stage)
    kept = (tmp_path / 'x.noncoding.fasta').read_text()
    assert kept == '>t2|lnc_RNA|S\nAAAA\n>t4|lnc_RNA|S\nGGGG\n', 'kept records verbatim, in order'
    ids = (tmp_path / 'x.coding_ids.txt').read_text().split()
    assert ids == coding_ids, 'input stage strips to the bare tid, projected keeps the id'
    record = json.loads((tmp_path / 'x.td2_coding.data.json').read_text())
    assert record['identity']['stage'] == stage
    assert record['data'] == {'n_in': 4, 'n_coding': 2, 'n_kept': 2, 'psauron_min': 0.5}, (
        'hits outside the FASTA ignored'
    )


def test_empty_pep_keeps_everything(tmp_path):
    fa, pep = write(tmp_path / 'in.fa', FASTA), write(tmp_path / 'p.pep', '')
    run(tmp_path, fa, pep)
    assert (tmp_path / 'x.noncoding.fasta').read_text() == FASTA


def test_malformed_score_warns_and_is_not_coding(tmp_path):
    pep = write(tmp_path / 'p.pep', '>t4|lnc_RNA|S.p1 G  ORF type:complete (+),psauron_score=n/a len:9\nM*\n')
    fa = write(tmp_path / 'in.fa', FASTA)
    r, _ = run(tmp_path, fa, pep)
    assert 'no parseable psauron_score for t4|lnc_RNA|S' in r.stderr
    assert (tmp_path / 'x.noncoding.fasta').read_text() == FASTA


def test_record(tmp_path):
    fa, pep = write(tmp_path / 'in.fa', FASTA), write(tmp_path / 'p.pep', PEP)
    args = ['--fasta', fa, '--pep', pep, '--prefix', tmp_path / 'x', '--stage', 'projected']
    identity = '{"species": "T", "side": "target", "target": "T", "track": "lnc_RNA"}'
    run_script('td2_coding_filter.py', *args, '--identity', identity, check=True)
    record = json.loads((tmp_path / 'x.td2_coding.data.json').read_text())
    assert record['section'] == 'td2_coding' and record['identity']['stage'] == 'projected'
    assert record['data'] == {'n_in': 4, 'n_coding': 2, 'n_kept': 2, 'psauron_min': 0.5}
