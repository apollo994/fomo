#!/usr/bin/env python3
"""Check that the documentation points at things that exist: every relative Markdown link, and
every repository path written in backticks, in CLAUDE.md, README.md, CONTRIBUTING.md, CHANGELOG.md,
docs/ and plans/README.md. Run from the repository root; exit 1 lists what is missing.

A backticked path counts when it starts with a top-level directory of the repository and holds
no placeholder (`<...>`, `{...}`, `*`, `$`).
"""

import re
import sys
from pathlib import Path

DOCS = ['CLAUDE.md', 'README.md', 'CONTRIBUTING.md', 'CHANGELOG.md', 'plans/README.md']
TOP = ('assets/', 'bin/', 'conf/', 'dev/', 'docs/', 'modules/', 'plans/', 'subworkflows/', 'tests/', 'workflows/')
LINK = re.compile(r'\[[^\]]*\]\(([^)\s]+)\)')
CODE = re.compile(r'`([^`\s]+)`')
PLACEHOLDER = re.compile(r'[<>{}*$]')


def doc_files(root: Path) -> list[Path]:
    files = [root / f for f in DOCS if (root / f).exists()]
    return files + sorted((root / 'docs').glob('*.md'))


def missing(root: Path, doc: Path) -> list[str]:
    out = []
    text = doc.read_text(encoding='utf-8')
    for target in LINK.findall(text):
        if re.match(r'[a-z]+:', target) or target.startswith('#'):
            continue
        path = target.split('#', 1)[0]
        if path and not (doc.parent / path).exists():
            out.append(f'{doc.relative_to(root)}: link {target}')
    for span in CODE.findall(text):
        path = span.rstrip('/.,:;')
        if path.startswith(TOP) and not PLACEHOLDER.search(path) and not (root / path).exists():
            out.append(f'{doc.relative_to(root)}: path {span}')
    return out


def main() -> int:
    root = Path.cwd()
    problems = [m for doc in doc_files(root) for m in missing(root, doc)]
    for problem in problems:
        print(problem, file=sys.stderr)
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
