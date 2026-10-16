"""dev/: the link check, the comment budget, the plans index and the end-of-file fixer."""

import importlib.util
from pathlib import Path

DEV = Path(__file__).resolve().parents[2] / 'dev'


def load(name):
    spec = importlib.util.spec_from_file_location(name, DEV / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


check_links = load('check_links')
code_metrics = load('code_metrics')
plans_index = load('plans_index')
end_of_file = load('end_of_file')


def test_links_and_paths(tmp_path):
    (tmp_path / 'bin').mkdir()
    (tmp_path / 'bin/here.py').write_text('')
    (tmp_path / 'docs').mkdir()
    doc = tmp_path / 'docs/a.md'
    doc.write_text(
        '[ok](../bin/here.py) [web](https://example.org) [anchor](#x) [gone](missing.md)\n'
        '`bin/here.py` `bin/gone.py` `bin/<name>.py` `not/a/top/path`\n'
    )
    assert check_links.missing(tmp_path, doc) == ['docs/a.md: link missing.md', 'docs/a.md: path bin/gone.py']


def test_budget_phrases():
    history, incident = code_metrics.BUDGET_HISTORY, code_metrics.BUDGET_INCIDENT
    assert history.search('this is no longer needed') and history.search('the old parser')
    assert not history.search('models whose transcripts were dropped')
    assert incident.search('as seen for real on run agitated_x') and not incident.search('verified against gffcompare')


def test_repository_meets_its_budget():
    assert code_metrics.check() == 0


def test_plans_index_row(tmp_path):
    plan = tmp_path / '99_x.md'
    plan.write_text('# 99 — Use `bin/x.py`\n\n> **Status:** Implemented · **Date:** 2026-01-02\n> **Decision:** y\n')
    assert plans_index.row(plan) == '| [99_x.md](99_x.md) | 99 — Use bin/x.py | Implemented | 2026-01-02 |'
    plan.write_text('# 27z — Review\n\n## Status: done (2026-10-15)\n')
    assert plans_index.row(plan) == '| [99_x.md](99_x.md) | 27z — Review | done (2026-10-15) |  |'


def test_end_of_file_fixer(tmp_path, capsys):
    cases = {'ok.txt': b'a\n', 'none.txt': b'a', 'extra.txt': b'a\n\n\r\n', 'blank.txt': b'\n\n', 'empty.txt': b''}
    for name, data in cases.items():
        (tmp_path / name).write_bytes(data)
    assert end_of_file.main([str(tmp_path / n) for n in cases]) == 1
    assert {n: (tmp_path / n).read_bytes() for n in cases} == {
        'ok.txt': b'a\n',
        'none.txt': b'a\n',
        'extra.txt': b'a\n',
        'blank.txt': b'',
        'empty.txt': b'',
    }
    assert capsys.readouterr().out.count('fixed end of file') == 3
    assert end_of_file.main([str(tmp_path / 'ok.txt')]) == 0
