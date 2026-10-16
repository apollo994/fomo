"""bin/gffcompare_batch.sh and bin/gffcompare_stats_record.sh: one gffcompare per manifest row,
published names kept, one accuracy record per run, equal to fomo_stats' parser."""

import json
import subprocess
from pathlib import Path

import fomo_stats
import pytest

from conftest import BIN, needs, run_script, write

FIXTURES = Path(__file__).resolve().parents[1] / 'records' / 'fixtures'
REF = (
    'c\tref\tgene\t100\t1000\t.\t+\t.\tID=g1\n'
    'c\tref\tlnc_RNA\t100\t1000\t.\t+\t.\tID=r1;Parent=g1\n'
    'c\tref\texon\t100\t200\t.\t+\t.\tParent=r1\n'
    'c\tref\texon\t500\t1000\t.\t+\t.\tParent=r1\n'
)
IDENTITY = '{"species":"T","side":"target","target":"T","track":"lnc_RNA"}'


def query(tid, src):
    uid = f'{tid}|{src}'
    return (
        f'c\tfomo\tgene\t100\t1000\t.\t+\t.\tID=gene-{uid};source={src}\n'
        f'c\tfomo\ttranscript\t100\t1000\t.\t+\t.\tID={uid};Parent=gene-{uid};source={src}\n'
        f'c\tfomo\texon\t100\t200\t.\t+\t.\tID={uid}.exon1;Parent={uid}\n'
        f'c\tfomo\texon\t500\t1000\t.\t+\t.\tID={uid}.exon2;Parent={uid}\n'
    )


def batch(tmp_path, *extra, path=None):
    write(tmp_path / 'ref.gff3', REF)
    write(tmp_path / 'T.allModels.lnc_RNA.raw.gff3', query('t1', 'A'))
    write(tmp_path / 'T.from_A.lnc_RNA.projected.gff3', query('t1', 'A'))
    write(
        tmp_path / 'm.tsv',
        'file\tsource\nT.allModels.lnc_RNA.raw.gff3\tallModels_raw\nT.from_A.lnc_RNA.projected.gff3\tA\n',
    )
    cmd = [str(BIN / 'gffcompare_batch.sh'), '--manifest', 'm.tsv', '--reference', 'ref.gff3', '--target', 'T']
    env = None if path is None else {'PATH': f'{path}:{BIN}:/usr/bin:/bin'}
    return subprocess.run(cmd + ['--track', 'lnc_RNA', *extra], cwd=tmp_path, capture_output=True, text=True, env=env)


@needs('gffcompare')
def test_outputs_and_records(tmp_path):
    r = batch(tmp_path, '--identity', IDENTITY, '--', '-M', '--no-merge')
    assert r.returncode == 0, r.stderr
    for source in ('allModels_raw', 'A'):
        prefix = f'T.from_{source}.lnc_RNA'
        for suffix in ('gffcompare.stats', 'gffcompare.tracking', 'gffcompare.loci'):
            assert (tmp_path / f'{prefix}.{suffix}').exists(), suffix
        assert (tmp_path / f'{prefix}.gffcompare.{prefix}.projected.gff3.tmap').exists(), 'tmap named after the query'
        record = json.loads((tmp_path / f'{prefix}.accuracy.data.json').read_text())
        assert record['identity'] == json.loads(IDENTITY) | {'source': source}
        assert record['data']['levels']['Transcript'] == {'sn': 100.0, 'pr': 100.0}


def test_failing_gffcompare_fails_the_batch(tmp_path):
    fake = tmp_path / 'fake'
    fake.mkdir()
    write(fake / 'gffcompare', '#!/bin/sh\nexit 3\n').chmod(0o755)
    assert batch(tmp_path, path=fake).returncode != 0


@pytest.mark.parametrize('stats', sorted(FIXTURES.glob('*.gffcompare.stats')), ids=lambda p: p.name)
def test_record_equals_the_python_parser(stats):
    r = run_script('gffcompare_stats_record.sh', IDENTITY, stats)
    assert r.returncode == 0, r.stderr
    record = json.loads(r.stdout)
    assert record['section'] == 'accuracy'
    assert record['data'] == fomo_stats.parse_gffcompare_stats(str(stats))
