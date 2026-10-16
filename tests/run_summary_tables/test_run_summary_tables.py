"""bin/run_summary_tables.py: regression against the summary tables of real test-profile runs (golden V1, V2).

The fixtures are the data records of those runs (every *.data.json their tasks wrote) and the
tables RUN_SUMMARY_TABLES published.
"""

import json
import shutil
from pathlib import Path

import pytest

from conftest import run_script

FIXTURES = Path(__file__).parent / 'fixtures'
ARGS = {
    'V1': ['--include-mrna', 'false', '--include-decoy', 'false', '--psauron-min', '0.5'],
    'V2': ['--include-mrna', 'true', '--include-decoy', 'true', '--psauron-min', '0.5'],
}


def run(fx, records, out, monkeypatch, variant):
    monkeypatch.chdir(out)
    r = run_script('run_summary_tables.py', '--records', records, '--roles', fx / 'species_roles.csv', *ARGS[variant])
    assert r.returncode == 0, r.stderr


@pytest.mark.parametrize('variant', sorted(ARGS))
def test_tables_match_golden_run(variant, tmp_path, monkeypatch):
    fx = FIXTURES / variant
    run(fx, fx / 'records', tmp_path, monkeypatch, variant)
    expected = sorted(p.name for p in (fx / 'expected').iterdir())
    assert sorted(p.name for p in tmp_path.iterdir()) == expected
    for name in expected:
        assert (tmp_path / name).read_bytes() == (fx / 'expected' / name).read_bytes(), name


def test_curated_genes_come_from_the_novel_table(tmp_path, monkeypatch):
    """observed and novel share metric names (e.g. `genes`): the curated count is the novel one."""
    fx = FIXTURES / 'V1'
    records = tmp_path / 'records'
    shutil.copytree(fx / 'records', records)
    for path in records.glob('*.curation.data.json'):
        record = json.loads(path.read_text())
        record['data']['observed']['genes'] = 999999
        path.write_text(json.dumps(record))
    out = tmp_path / 'out'
    out.mkdir()
    run(fx, records, out, monkeypatch, 'V1')
    assert (out / 'run_target_types_mqc.tsv').read_bytes() == (fx / 'expected/run_target_types_mqc.tsv').read_bytes()
