"""Repository conventions that no single script test covers: module shape, container labels and
the vendored binary."""

import hashlib
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
MODULES = sorted((ROOT / 'modules/local').glob('*.nf'))
CONTAINERS = (ROOT / 'conf/containers.config').read_text()

# MULTIQC consumes the versions topic, so it cannot emit to it.
NO_VERSIONS = {'multiqc.nf'}
# Gated by `ext.when` in conf/modules.config.
WHEN_ALLOWED = {'curate_models.nf'}


def labels(text):
    return re.findall(r"^\s*label\s+'([^']+)'", text, re.M)


@pytest.mark.parametrize('module', MODULES, ids=lambda p: p.name)
def test_module_shape(module):
    text = module.read_text()
    process = re.search(r'^process\s+(\w+)', text, re.M).group(1)
    assert module.stem == process.lower(), 'the file name is the process name in lower case'
    img = [lbl for lbl in labels(text) if lbl.startswith('img_')]
    assert len(img) == 1, f'exactly one img_* label, found {img}'
    assert not re.search(r'^\s*(container|conda)\b', text, re.M), 'images live in conf/containers.config'
    assert re.search(r'^\s*stub:', text, re.M), 'every module has a stub'
    assert ('when:' in text) == (module.name in WHEN_ALLOWED)
    assert ('topic: versions' in text) == (module.name not in NO_VERSIONS)


def test_img_labels_match_containers_config():
    used = {lbl for m in MODULES for lbl in labels(m.read_text()) if lbl.startswith('img_')}
    defined = set(re.findall(r"withLabel:\s*'(img_\w+)'", CONTAINERS))
    assert used == defined


def test_gff_feature_stats_modules_are_amd64_only():
    for m in MODULES:
        text = m.read_text()
        runs_it = re.search(r'^\s+gff-feature-stats ', text, re.M) is not None
        assert runs_it == ('amd64_only' in labels(text)), m.name


def test_no_nf_core_modules():
    assert not (ROOT / 'modules/nf-core').exists()


def test_gff_feature_stats_checksum():
    expected = (ROOT / 'bin/gff-feature-stats.sha256').read_text().split()[0]
    assert hashlib.sha256((ROOT / 'bin/gff-feature-stats').read_bytes()).hexdigest() == expected
