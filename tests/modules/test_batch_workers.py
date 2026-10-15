"""The per-target batch modules run one worker per file through xargs: a failing command inside
a worker must fail the worker (strict mode), not be masked by the commands after it."""

import re
import subprocess
from pathlib import Path

import pytest

MODULES = Path(__file__).resolve().parents[2] / 'modules/local'

CASES = {
    'gffcompare_batch.nf': 'gffcompare() { return 1; }\nmv() { return 0; }\n',
    'gff_stats_projected_batch.nf': 'gff-feature-stats() { return 1; }\ngff_stats_to_mqc.py() { return 0; }\n',
}


def worker(module):
    """The bash run_one() of a module's script block, with Groovy interpolations neutralised."""
    text = (MODULES / module).read_text()
    body = re.search(r'\n(\s*run_one\(\) \{\n.*?\n\s*\}\n)', text, re.S).group(1)
    body = re.sub(r'(?<!\\)\$\{[^}]*\}', 'X', body)  # Groovy ${...}
    return body.replace('\\$', '$').replace('\\\\', '\\')


@pytest.mark.parametrize('module', sorted(CASES))
def test_failing_first_command_fails_the_worker(module, tmp_path):
    script = CASES[module] + worker(module) + '\nrun_one T.from_S.lnc_RNA.projected.gff3\n'
    r = subprocess.run(['bash', '-c', script], cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode != 0, 'a failing gffcompare / gff-feature-stats must fail the worker'
