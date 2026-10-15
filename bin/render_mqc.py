#!/usr/bin/env python3
"""
Render the custom content of one target's MultiQC report from data records.

Every section of fomo_sections.SECTIONS becomes one TSV per record (matched by its `sp` pattern
in assets/multiqc/sections.yml), and the accuracy records become one scatter JSON. A record is the
target's when its identity names that target, and everyone's (source-side) when it names none.
Runs inside the MULTIQC task, before multiqc; stdlib only.

Usage: render_mqc.py --records DIR --target T --out DIR
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fomo_io  # noqa: E402
import fomo_sections  # noqa: E402
import fomo_stats  # noqa: E402

log = fomo_io.get_logger('render_mqc')


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--records', required=True, type=Path, help='directory of *.data.json records')
    p.add_argument('--target', required=True, help='the target whose report this is')
    p.add_argument('--out', required=True, type=Path, help='output directory')
    return p.parse_args()


def write_table(path: Path, rows: list[dict]) -> None:
    columns = list(rows[0])
    with open(path, 'w', encoding='utf-8') as out:
        out.write('\t'.join(columns) + '\n')
        for row in rows:
            out.write('\t'.join(row[c] for c in columns) + '\n')


def render(records: list[dict], target: str, out: Path) -> int:
    """Write every section's tables and the scatter; -> files written."""
    mine = [r for r in records if r['identity'].get('target') in (None, target)]
    n_files = 0
    for section in fomo_sections.SECTIONS:
        for rec in mine:
            if rec['section'] == section.record and section.select(rec['identity']):
                stem, rows = section.rows(rec)
                write_table(out / f'{stem}{section.suffix}', rows)
                n_files += 1
    accuracy = [r for r in mine if r['section'] == 'accuracy']
    if accuracy:
        with open(out / f'{target}{fomo_sections.SCATTER_SUFFIX}', 'w', encoding='utf-8') as fh:
            json.dump(fomo_sections.scatter(accuracy), fh, indent=2)
        n_files += 1
    return n_files


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    n_files = render(fomo_stats.read_records(args.records), args.target, args.out)
    log.info(f'{args.target}: {n_files} MultiQC files')
    return 0


if __name__ == '__main__':
    fomo_io.run(main)
