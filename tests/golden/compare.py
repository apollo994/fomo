#!/usr/bin/env python3
"""Compare two deliverable snapshots (snapshot.py). Exits 1 on any missing, extra or changed path.

Usage: compare.py EXPECTED.tsv ACTUAL.tsv
"""

import sys


def load(path: str) -> dict:
    with open(path, encoding='utf-8') as fh:
        next(fh)
        return {p: (rule, sha) for p, rule, sha in (line.rstrip('\n').split('\t') for line in fh)}


def main() -> int:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    expected, actual = load(sys.argv[1]), load(sys.argv[2])
    missing = sorted(set(expected) - set(actual))
    extra = sorted(set(actual) - set(expected))
    changed = sorted(p for p in set(expected) & set(actual) if expected[p] != actual[p])
    for label, paths in (('missing', missing), ('extra', extra), ('changed', changed)):
        for p in paths:
            print(f'{label}\t{p}')
    n = len(missing) + len(extra) + len(changed)
    print(
        f'# {len(expected)} expected, {len(actual)} actual: {len(missing)} missing, {len(extra)} extra, '
        f'{len(changed)} changed',
        file=sys.stderr,
    )
    return 1 if n else 0


if __name__ == '__main__':
    sys.exit(main())
