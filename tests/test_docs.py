"""docs/data_records.md documents exactly the record sections of fomo_stats.SECTIONS, each with its
required identity keys."""

import re
from pathlib import Path

import fomo_stats

DOC = (Path(__file__).resolve().parent.parent / 'docs/data_records.md').read_text()


def test_every_section_and_only_those():
    assert re.findall(r'^### (\w+)$', DOC, re.M) == list(fomo_stats.SECTIONS)


def test_required_identity_keys():
    for section, keys in fomo_stats.SECTIONS.items():
        body = DOC.split(f'### {section}\n', 1)[1].split('\n### ', 1)[0]
        required = re.search(r'Requires\s+(.*?)\.', body, re.S).group(1)
        assert re.findall(r'`(\w+)`', required) == list(keys), section
