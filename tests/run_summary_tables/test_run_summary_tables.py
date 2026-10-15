"""bin/run_summary_tables.py: regression against the summary tables of real test-profile runs (golden V1, V2).

The fixtures are the staged inputs of the RUN_SUMMARY_TABLES task of those runs (pruned to the files
the script reads) and the tables it published.
"""

from pathlib import Path

import pytest

from conftest import run_script

FIXTURES = Path(__file__).parent / 'fixtures'
ARGS = {
    'V1': ['--include-mrna', 'false', '--include-decoy', 'false', '--psauron-min', '0.5'],
    'V2': ['--include-mrna', 'true', '--include-decoy', 'true', '--psauron-min', '0.5'],
}


@pytest.mark.parametrize('variant', sorted(ARGS))
def test_tables_match_golden_run(variant, tmp_path, monkeypatch):
    fx = FIXTURES / variant
    monkeypatch.chdir(tmp_path)
    r = run_script(
        'run_summary_tables.py', '--mqc-dir', fx / 'mqc', '--roles', fx / 'species_roles.csv', *ARGS[variant]
    )
    assert r.returncode == 0, r.stderr
    expected = sorted(p.name for p in (fx / 'expected').iterdir())
    assert sorted(p.name for p in tmp_path.iterdir()) == expected
    for name in expected:
        assert (tmp_path / name).read_bytes() == (fx / 'expected' / name).read_bytes(), name
