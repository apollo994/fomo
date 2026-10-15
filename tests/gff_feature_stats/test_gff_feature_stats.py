"""The vendored bin/gff-feature-stats binary runs and reads a gzipped GFF3."""

import json
import platform
import subprocess

import pytest

from conftest import BIN

pytestmark = pytest.mark.skipif(
    platform.system() != 'Linux' or platform.machine() != 'x86_64', reason='x86-64 Linux binary'
)
BINARY = BIN / 'gff-feature-stats'


def test_version():
    r = subprocess.run([BINARY, '-V'], capture_output=True, text=True)
    assert r.returncode == 0 and r.stdout.strip().startswith('gff-feature-stats 0.2.0')


def test_reads_gzipped_annotation(tmp_path):
    gff = BIN.parent / 'assets/test_data/Lycaena_thersamon_265360/GCA_965194355.1/Lycaena_thersamon_265360.gff3.gz'
    out = tmp_path / 'stats.json'
    subprocess.run([BINARY, gff, out], check=True, capture_output=True)
    assert json.loads(out.read_text())
