"""bin/samtools_stats_record.sh and bin/seqkit_stats_record.sh: records carry the values of the
alignment and SeqKit tables, read from captured tool outputs of a V2 run."""

import json
from pathlib import Path

import pytest

from conftest import run_script

FIXTURES = Path(__file__).parent / 'fixtures' / 'emitters'
IDENTITY = '{"target": "T", "side": "target", "track": "lnc_RNA"}'


def record(script, path):
    r = run_script(script, IDENTITY, path)
    assert r.returncode == 0, r.stderr
    return json.loads(r.stdout)


@pytest.mark.parametrize('stats', sorted(FIXTURES.glob('*.stats')), ids=lambda p: p.name)
def test_alignment_record_matches_the_adapter(stats, tmp_path):
    rec = record('samtools_stats_record.sh', stats)
    assert rec['schema'] == 'fomo-data/1' and rec['section'] == 'alignment'
    assert rec['identity'] == json.loads(IDENTITY)
    out = tmp_path / 'a_samtools_align_mqc.tsv'
    run_script('samtools_to_mqc.py', '--input', stats, '--name', 'x', '--output', out, check=True)
    header, row = (line.split('\t') for line in out.read_text().splitlines())
    for key, value in zip(header[1:], row[1:], strict=True):
        assert rec['data'][key] == pytest.approx(float(value)), key


@pytest.mark.parametrize('tsv', sorted(FIXTURES.glob('*.tsv')), ids=lambda p: p.name)
def test_fasta_stats_record_keeps_every_column(tsv):
    rec = record('seqkit_stats_record.sh', tsv)
    assert rec['section'] == 'fasta_stats'
    header, row = (line.split('\t') for line in tsv.read_text().splitlines())
    for key, value in zip(header, row, strict=True):
        if key == 'file':
            assert key not in rec['data']
        elif key in ('format', 'type'):
            assert rec['data'][key] == value
        else:
            assert rec['data'][key] == pytest.approx(float(value)), key


def test_empty_inputs_still_give_a_record(tmp_path):
    empty = tmp_path / 'empty'
    empty.write_text('')
    assert record('samtools_stats_record.sh', empty)['data']['reads_mapped'] is None
    assert record('seqkit_stats_record.sh', empty)['data'] == {}
