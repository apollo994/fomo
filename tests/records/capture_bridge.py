#!/usr/bin/env python3
"""Capture the bridge-test fixture from a finished run: every MultiQC table produced from a
measurement, paired with the data record of the same measurement through task lineage (an
adapter's staged input links back to the task that produced it), never through Sample labels.

Usage: capture_bridge.py <run work dir> <out.json>
"""

import json
import sys
from pathlib import Path

ADAPTER_SUFFIXES = {
    '_gffstats_input_mqc.tsv': 'gffstats_input',
    '_gffstats_genes_mqc.tsv': 'gffstats_genes',
    '_samtools_align_mqc.tsv': 'alignment',
    '_seqkit_mqc.tsv': 'fasta_stats',
}
SAME_TASK = {
    '_transcript_filter_mqc.tsv': ('transcript_filter', '.transcript_filter.data.json'),
    '_td2_coding_mqc.tsv': ('td2_coding', '.td2_coding.data.json'),
    '_gffstats_projection_mqc.tsv': ('gffstats_projection', '.gffstats.data.json'),
}


def load(path):
    return json.loads(path.read_text())


def record_of(task_dir, suffix):
    found = sorted(task_dir.glob(f'*{suffix}'))
    if len(found) != 1:
        raise SystemExit(f'{task_dir}: expected one *{suffix}, found {len(found)}')
    return load(found[0])


def adapter_pairs(task, tsv, kind):
    """An adapter task: its staged input is a link into the producing task."""
    inputs = [p for p in task.iterdir() if p.is_symlink() and not p.name.startswith('.')]
    if len(inputs) != 1:
        raise SystemExit(f'{task}: expected one staged input, found {len(inputs)}')
    producer = inputs[0].resolve().parent
    suffix = {'alignment': '.alignment.data.json', 'fasta_stats': '.fasta_stats.data.json'}.get(
        kind, '.gffstats.data.json'
    )
    return {'kind': kind, 'record': record_of(producer, suffix), 'tables': {tsv.name: tsv.read_text()}}


def produced(task, pattern):
    """Files the task wrote itself (staged inputs are symlinks)."""
    return sorted(p for p in task.glob(pattern) if not p.is_symlink())


def task_pairs(task):
    pairs = []
    for tsv in produced(task, '*_mqc.tsv'):
        kind = next((k for s, k in ADAPTER_SUFFIXES.items() if tsv.name.endswith(s)), None)
        if kind:
            pairs.append(adapter_pairs(task, tsv, kind))
            continue
        for suffix, (same_kind, record_suffix) in SAME_TASK.items():
            if tsv.name.endswith(suffix):
                record = task / (tsv.name[: -len(suffix)] + record_suffix)
                pairs.append({'kind': same_kind, 'record': load(record), 'tables': {tsv.name: tsv.read_text()}})
    curation = produced(task, '*_curation_*_mqc.tsv')
    if curation:
        record = record_of(task, '.curation.data.json')
        pairs = [p for p in pairs if not p['tables'].keys() & {t.name for t in curation}]
        pairs.append({'kind': 'curation', 'record': record, 'tables': {t.name: t.read_text() for t in curation}})
    for stats in produced(task, '*.gffcompare.stats'):
        record = task / stats.name.replace('.gffcompare.stats', '.accuracy.data.json')
        pairs.append({'kind': 'accuracy', 'record': load(record), 'tables': {stats.name: stats.read_text()}})
    return pairs


def main():
    work, out = Path(sys.argv[1]), Path(sys.argv[2])
    pairs = [p for task in sorted(work.glob('*/*')) if task.is_dir() for p in task_pairs(task)]
    out.write_text(json.dumps(pairs, indent=1) + '\n')
    kinds = {}
    for p in pairs:
        kinds[p['kind']] = kinds.get(p['kind'], 0) + 1
    print(f'{len(pairs)} pairs: {kinds}')


if __name__ == '__main__':
    main()
