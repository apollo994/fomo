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


def test_curated_genes_come_from_the_novel_table(tmp_path, monkeypatch):
    """observed and novel curation rows share column names (e.g. `genes`): the curated count is
    the novel one, whatever order the tables are read in."""
    fx = FIXTURES / 'V1'
    mqc = tmp_path / 'mqc'
    mqc.mkdir()
    for f in (fx / 'mqc').iterdir():
        text = f.read_text()
        if f.name.endswith('_curation_observed_mqc.tsv'):
            head, row = text.splitlines()[:2]
            cols, vals = head.split('\t'), row.split('\t')
            vals[cols.index('genes')] = '999999'
            text = head + '\n' + '\t'.join(vals) + '\n'
        (mqc / f.name).write_text(text)
    out = tmp_path / 'out'
    out.mkdir()
    monkeypatch.chdir(out)
    r = run_script('run_summary_tables.py', '--mqc-dir', mqc, '--roles', fx / 'species_roles.csv', *ARGS['V1'])
    assert r.returncode == 0, r.stderr
    assert (out / 'run_target_types_mqc.tsv').read_bytes() == (fx / 'expected/run_target_types_mqc.tsv').read_bytes()
