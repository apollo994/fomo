"""bin/samtools_stats_record.sh and bin/seqkit_stats_record.sh: records carry the tool values, read
from captured tool outputs of a V2 run."""

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


SN_KEYS = {
    'raw_total_sequences': 'raw total sequences',
    'reads_mapped': 'reads mapped',
    'non_primary_alignments': 'non-primary alignments',
    'error_rate': 'error rate',
    'average_quality': 'average quality',
    'average_length': 'average length',
}


@pytest.mark.parametrize('stats', sorted(FIXTURES.glob('*.stats')), ids=lambda p: p.name)
def test_alignment_record_carries_the_sn_fields(stats):
    rec = record('samtools_stats_record.sh', stats)
    assert rec['schema'] == 'fomo-data/1' and rec['section'] == 'alignment'
    assert rec['identity'] == json.loads(IDENTITY)
    sn = {}
    for line in stats.read_text().splitlines():
        if line.startswith('SN\t'):
            _, key, value = line.split('\t')[:3]
            sn[key.rstrip(':')] = float(value.split('#')[0])
    for key, sn_key in SN_KEYS.items():
        assert rec['data'][key] == pytest.approx(sn[sn_key]), key
    mapped = 100 * sn['reads mapped'] / sn['raw total sequences']
    assert rec['data']['reads_mapped_percent'] == pytest.approx(round(mapped, 2))


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
