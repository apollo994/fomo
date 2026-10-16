"""bin/render_mqc.py and bin/fomo_sections.py: the rendered MultiQC content equals the tables the
adapter processes wrote, row label for row label and value for value.

Fixtures: bridge_V{1,2}.json, the adapter tables of the V1 and V2 golden runs paired with the
record of the same measurement through task lineage (captured before the adapters were deleted).
"""

import csv
import io
import json
from pathlib import Path

import fomo_sections
import pytest
import render_mqc

FIXTURES = Path(__file__).parent / 'fixtures'
PAIRS = [
    pytest.param(pair, id=f'{variant}:{pair["kind"]}:{next(iter(pair["tables"]))}')
    for variant in ('V1', 'V2')
    for pair in json.loads((FIXTURES / f'bridge_{variant}.json').read_text())
]
SECTION_OF_KIND = {
    'gffstats_input': 'gffstats_input',
    'gffstats_genes': 'gffstats_genes',
    'gffstats_projection': 'gffstats_projection',
    'fasta_stats': 'seqkit_stats',
    'alignment': 'samtools_align',
    'transcript_filter': 'transcript_filter',
    'td2_coding': 'td2_coding',
}
SECTION = {s.id: s for s in fomo_sections.SECTIONS}


def rows(text):
    return list(csv.DictReader(io.StringIO(text), delimiter='\t'))


def same(old: str, new: str) -> bool:
    if old == new:
        return True
    try:
        return float(old) == pytest.approx(float(new))
    except ValueError:
        return False


def assert_same_rows(old_rows, new_rows):
    assert [r['Sample'] for r in new_rows] == [r['Sample'] for r in old_rows]
    for old, new in zip(old_rows, new_rows, strict=True):
        assert list(new) == list(old), 'same columns, same order'
        for key in old:
            assert same(old[key], new[key]), (old['Sample'], key, old[key], new[key])


def rendered(section_id, record):
    section = SECTION[section_id]
    assert section.select(record['identity'])
    return section.rows(record)[1]


@pytest.mark.parametrize('pair', PAIRS)
def test_rendered_content_equals_the_adapter_table(pair):
    kind, record = pair['kind'], pair['record']
    if kind in SECTION_OF_KIND:
        (text,) = pair['tables'].values()
        assert_same_rows(rows(text), rendered(SECTION_OF_KIND[kind], record))
    elif kind == 'curation':
        for name, text in pair['tables'].items():
            table = name.rsplit('_curation_', 1)[1].removesuffix('_mqc.tsv')
            assert_same_rows(rows(text), rendered(f'curation_{table}', record))
    elif kind == 'accuracy':
        (text,) = pair['tables'].values()
        levels = {}
        for line in text.splitlines():
            if 'level:' in line:
                toks = line.replace('|', '').replace('Intron chain', 'Intron_chain').split()
                levels[toks[0]] = (round(float(toks[2]) / 100, 4), round(float(toks[3]) / 100, 4))
        plot = fomo_sections.scatter([record])
        name = f'{record["identity"]["source"]}.{record["identity"]["track"]}'
        assert [(lvl[name]['x'], lvl[name]['y']) for lvl in plot['data']] == list(levels.values())
    else:
        raise AssertionError(kind)


def test_render_writes_the_target_and_broadcast_content(tmp_path):
    unique = {json.dumps(p.values[0]['record'], sort_keys=True) for p in PAIRS if p.id.startswith('V1:')}
    records = [json.loads(r) for r in sorted(unique)]
    out = tmp_path / 'out'
    out.mkdir()
    n = render_mqc.render(records, 'Lycaena_hippothoe_580924', out)
    names = sorted(p.name for p in out.iterdir())
    assert n == len(names) and not any(name.startswith('Lycaena_alciphron_282377.from_') for name in names)
    assert any(name.endswith('_curation_novel_mqc.tsv') for name in names)
    assert any(name.startswith('Lycaena_thersamon_265360.') for name in names), 'source-side records are broadcast'


def test_colours_are_distinct_for_many_sources():
    colours = fomo_sections.source_colours({f'sp{i}' for i in range(40)} | {'allModels_raw'})
    assert len(set(colours.values())) == 41
