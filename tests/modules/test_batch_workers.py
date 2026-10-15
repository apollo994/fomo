"""GFF_STATS_PROJECTED_BATCH runs its script under strict mode over the manifest: a failing
gff-feature-stats fails the task, and each row's record carries that row's source."""

import json
import re
import subprocess
from pathlib import Path

MODULE = Path(__file__).resolve().parents[2] / 'modules/local/gff_stats_projected_batch.nf'
IDENTITY = '{"species":"T","side":"target","target":"T","track":"lnc_RNA","kind":"projected"}'


def script(tmp_path):
    """The module's script block as bash, Groovy interpolations replaced."""
    body = re.search(r'\n    script:\n.*?"""\n(.*?)\n    """', MODULE.read_text(), re.S).group(1)
    body = body.replace('${task.ext.identity}', IDENTITY).replace('${manifest}', 'manifest.tsv')
    body = re.sub(r'(?<!\\)\$\{[^}]*\}', 'X', body)
    (tmp_path / 'manifest.tsv').write_text('file\tsource\nT.allModels.lnc_RNA.raw.gff3\tallModels_raw\nb.gff3\tB\n')
    return body.replace('\\$', '$').replace('\\\\', '\\')


def run(tmp_path, stubs):
    return subprocess.run(['bash', '-c', stubs + script(tmp_path)], cwd=tmp_path, capture_output=True, text=True)


def test_failing_gff_feature_stats_fails_the_task(tmp_path):
    r = run(tmp_path, 'gff-feature-stats() { return 1; }\ngff_stats_to_mqc.py() { :; }\n')
    assert r.returncode != 0


def test_records_take_their_source_from_the_manifest(tmp_path):
    stubs = 'gff-feature-stats() { echo \'{"transcript_type_stats":{}}\' > "$2"; }\ngff_stats_to_mqc.py() { :; }\n'
    r = run(tmp_path, stubs)
    assert r.returncode == 0, r.stderr
    for source in ('allModels_raw', 'B'):
        record = json.loads((tmp_path / f'X.from_{source}.X.projected.gffstats.data.json').read_text())
        assert record['identity'] == json.loads(IDENTITY) | {'source': source}
        assert record['section'] == 'gffstats' and record['data'] == {'transcript_type_stats': {}}
