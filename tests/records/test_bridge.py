"""Bridge between the MultiQC tables and the data records: for every table a run produced, the
record of the same measurement (paired by task lineage in capture_bridge.py, not by Sample label)
holds the same values. This is what lets the table path be deleted.

Fixtures: capture_bridge.py on the work dirs of the V1 and V2 golden runs.
"""

import csv
import io
import json
from pathlib import Path

import curation_to_mqc
import fomo_stats
import gff_stats_to_mqc
import pytest

FIXTURES = Path(__file__).parent / 'fixtures'
PAIRS = [
    pytest.param(pair, id=f'{variant}:{pair["kind"]}:{next(iter(pair["tables"]))}')
    for variant in ('V1', 'V2')
    for pair in json.loads((FIXTURES / f'bridge_{variant}.json').read_text())
]


def rows(text):
    return list(csv.DictReader(io.StringIO(text), delimiter='\t'))


def only_table(pair):
    (text,) = pair['tables'].values()
    return text


def without_sample(row):
    return {k: v for k, v in row.items() if k != 'Sample'}


def numeric_equal(table_value, record_value):
    if table_value == 'NA':
        return record_value is None
    return record_value == pytest.approx(float(table_value))


def check_gffstats(pair, builder, role, kind):
    expected = builder(pair['record']['data'], role, kind)
    assert [without_sample(r) for r in rows(only_table(pair))] == expected


def check(pair):
    record, identity, data = pair['record'], pair['record']['identity'], pair['record']['data']
    fomo_stats.make_record(record['section'], identity, data)  # validates
    kind = pair['kind']
    if kind == 'gffstats_input':
        check_gffstats(pair, gff_stats_to_mqc.transcript_rows, identity['side'], identity['kind'])
    elif kind == 'gffstats_genes':
        check_gffstats(pair, gff_stats_to_mqc.gene_rows, identity['side'], identity['kind'])
    elif kind == 'gffstats_projection':
        check_gffstats(pair, gff_stats_to_mqc.transcript_rows, 'projection', 'projected')
        assert identity['kind'] == 'projected' and identity['source']
    elif kind in ('alignment', 'fasta_stats', 'transcript_filter'):
        (row,) = rows(only_table(pair))
        for key, value in without_sample(row).items():
            if key in ('format', 'type', 'feature'):
                continue
            assert numeric_equal(value, data[key]), key
    elif kind == 'td2_coding':
        (row,) = rows(only_table(pair))
        assert row['stage'] == identity['stage']
        assert [int(row[k]) for k in ('n_in', 'n_coding', 'n_kept')] == [
            data[k] for k in ('n_in', 'n_coding', 'n_kept')
        ]
    elif kind == 'curation':
        check_curation(pair, data)
    elif kind == 'accuracy':
        check_accuracy(pair, identity, data)
    else:
        raise AssertionError(kind)


def check_curation(pair, data):
    tables = {name.rsplit('_curation_', 1)[1].split('_')[0]: rows(text)[0] for name, text in pair['tables'].items()}
    for key in curation_to_mqc.OBSERVED:
        assert numeric_equal(tables['observed'][key], data['observed'].get(key, 0)), key
    assert tables['discarded']['reference'] == ('yes' if data['has_reference'] else 'no')
    for key in curation_to_mqc.DISCARDED:
        expected = data['discarded_by_reference'].get(key, 0) if data['has_reference'] else None
        assert numeric_equal(tables['discarded'][key], expected), key
    for key in curation_to_mqc.NOVEL:
        assert numeric_equal(tables['novel'][key], data['novel'].get(key, 0)), key
    for loc in curation_to_mqc.LOCATIONS:
        assert numeric_equal(tables['novel'][f'location:{loc}'], data['novel'].get(f'genes_location:{loc}', 0))


def check_accuracy(pair, identity, data):
    (name,) = pair['tables']
    parsed = fomo_stats.parse_identity(name)
    assert (parsed.target, parsed.source) == (identity['target'], identity['source'])
    assert parsed.feature_type + ('.decoy' if parsed.decoy else '') == identity['track']
    levels = {}
    for line in only_table(pair).splitlines():
        if 'level:' in line:
            toks = line.replace('|', '').replace('Intron chain', 'Intron_chain').split()
            levels[toks[0]] = {'sn': float(toks[2]), 'pr': float(toks[3])}
    assert data['levels'] == levels


@pytest.mark.parametrize('pair', PAIRS)
def test_record_carries_the_table(pair):
    check(pair)


def test_every_producer_is_covered():
    kinds = {pair.values[0]['kind'] for pair in PAIRS}
    assert kinds == {
        'gffstats_input',
        'gffstats_genes',
        'gffstats_projection',
        'alignment',
        'fasta_stats',
        'transcript_filter',
        'td2_coding',
        'curation',
        'accuracy',
    }
