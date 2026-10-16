"""GFF_STATS runs its script under strict mode over its rows (a manifest, or the one GFF3): a
failing gff-feature-stats fails the task, and each record carries the row's source."""

import json
import re
import subprocess
from pathlib import Path

MODULE = Path(__file__).resolve().parents[2] / 'modules/local/gff_stats.nf'
IDENTITY = '{"species":"T","side":"target","target":"T","track":"lnc_RNA","kind":"projected"}'
STATS = 'gff-feature-stats() { echo \'{"transcript_type_stats":{}}\' > "$2"; }\n'


def script(rows):
    """The module's script block as bash, Groovy interpolations replaced."""
    body = re.search(r'\n    script:\n.*?"""\n(.*?)\n    """', MODULE.read_text(), re.S).group(1)
    for groovy, value in {
        '${task.ext.identity}': IDENTITY,
        '${rows}': rows,
        '${prefix}': 'P',
        '${meta.target_id}': 'T',
        '${meta.track}': 'lnc_RNA',
    }.items():
        body = body.replace(groovy, value)
    assert not re.search(r'(?<!\\)\$\{[a-z]', body), 'an interpolation was not replaced'
    return body.replace('\\$', '$').replace('\\\\', '\\')


def run(tmp_path, stubs, rows):
    return subprocess.run(['bash', '-c', stubs + script(rows)], cwd=tmp_path, capture_output=True, text=True)


def manifest(tmp_path):
    (tmp_path / 'manifest.tsv').write_text('file\tsource\nT.allModels.lnc_RNA.raw.gff3\tallModels_raw\nb.gff3\tB\n')
    return 'tail -n +2 manifest.tsv'


def test_failing_gff_feature_stats_fails_the_task(tmp_path):
    assert run(tmp_path, 'gff-feature-stats() { return 1; }\n', manifest(tmp_path)).returncode != 0


def test_batch_records_take_their_source_from_the_manifest(tmp_path):
    r = run(tmp_path, STATS, manifest(tmp_path))
    assert r.returncode == 0, r.stderr
    for source in ('allModels_raw', 'B'):
        record = json.loads((tmp_path / f'T.from_{source}.lnc_RNA.projected.gffstats.data.json').read_text())
        assert record['identity'] == json.loads(IDENTITY) | {'source': source}
        assert record['section'] == 'gffstats' and record['data'] == {'transcript_type_stats': {}}


def test_single_file_record_takes_the_task_identity(tmp_path):
    r = run(tmp_path, STATS, "printf '%s\\t\\n' in.gff3")
    assert r.returncode == 0, r.stderr
    record = json.loads((tmp_path / 'P.gffstats.data.json').read_text())
    assert record['identity'] == json.loads(IDENTITY)
