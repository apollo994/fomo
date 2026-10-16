#!/usr/bin/env python3
"""Write plans/README.md: one row per plan with its title, status and date, read from the plan's
decision-record header (`> **Status:** ... · **Date:** ...`) or its `## Status:` line.
Run from the repository root after a plan's status changes."""

import re
from pathlib import Path

INTRO = """# Plans

Design decisions, one file per plan. Plans 00–26 carry a decision-record header (status, date,
decision); their bodies are kept as written. Plan 27 is the refactor that produced the current code
base: a master plan with subplans 27a–27k, a frozen copy of the plan as first approved in
[27_refactor_original/](27_refactor_original/), and the review of every deviation in
`27z_deviation_review.md` once 27k is done.

| plan | title | status | date |
|:---|:---|:---|:---|
"""


def row(path: Path) -> str:
    text = path.read_text(encoding='utf-8')
    title = text.split('\n', 1)[0].lstrip('# ').strip().replace('`', '')
    adr = re.search(r'^> \*\*Status:\*\* (.*?) · \*\*Date:\*\* (\S+)', text, re.M)
    if adr:
        status, date = adr.groups()
    else:
        line = re.search(r'^## Status: (.*)$', text, re.M)
        status, date = (line.group(1) if line else ''), ''
    return f'| [{path.name}]({path.name}) | {title} | {status} | {date} |'


def main() -> None:
    plans = sorted(p for p in Path('plans').glob('*.md') if p.name != 'README.md')
    Path('plans/README.md').write_text(INTRO + '\n'.join(row(p) for p in plans) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
