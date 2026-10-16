"""The MultiQC configuration agrees with what the pipeline writes: every rendered section has an
`sp` pattern and a place in the report order, every pattern belongs to a section, and the same
for the run summary against run_summary_tables.py's outputs."""

import fnmatch
import json
import shutil
import subprocess
import sys
from pathlib import Path

import fomo_sections
import yaml

from conftest import BIN

ASSETS = Path(__file__).resolve().parent.parent / 'assets/multiqc'
RUN_FIXTURE = Path(__file__).resolve().parent / 'run_summary_tables/fixtures/V2'


def load(name):
    return yaml.safe_load((ASSETS / name).read_text())


def test_per_target_sections_match_the_registry():
    sections, main = load('sections.yml'), load('main.yml')
    registry = {s.id: s for s in fomo_sections.SECTIONS}
    assert set(sections['sp']) == set(registry) == set(sections['custom_data'])
    for sid, section in registry.items():
        assert sections['sp'][sid]['fn'] == f'*{section.suffix}', sid
    ordered = set(main['report_section_order'])
    assert set(registry) | {fomo_sections.SCATTER_ID} <= ordered


def test_run_summary_sections_match_its_outputs(tmp_path):
    shutil.copytree(RUN_FIXTURE / 'records', tmp_path / 'records')
    args = ['--records', 'records', '--roles', str(RUN_FIXTURE / 'species_roles.csv'), '--include-mrna', 'true']
    subprocess.run(
        [sys.executable, str(BIN / 'run_summary_tables.py'), *args, '--include-decoy', 'true'],
        cwd=tmp_path,
        check=True,
        capture_output=True,
    )
    sections, main = load('run_summary_sections.yml'), load('run_summary_main.yml')
    tsvs = sorted(p.name for p in tmp_path.glob('*_mqc.tsv'))
    patterns = {sid: sp['fn'] for sid, sp in sections['sp'].items()}
    for name in tsvs:
        assert sum(fnmatch.fnmatch(name, fn) for fn in patterns.values()) == 1, name
    for sid, fn in patterns.items():
        assert any(fnmatch.fnmatch(name, fn) for name in tsvs), sid
    plot_ids = {json.loads(p.read_text())['id'] for p in tmp_path.glob('*_mqc.json')}
    assert set(patterns) | plot_ids == set(main['report_section_order'])
