#!/usr/bin/env python3
"""End every text file with exactly one newline (pre-commit hook). Files that need it are fixed in
place and reported, and the exit status is 1, as pre-commit's own end-of-file-fixer does.

Usage: end_of_file.py FILE...
"""

import sys


def fix(data: bytes) -> bytes:
    stripped = data.rstrip(b'\r\n')
    return stripped + b'\n' if stripped else b''


def main(paths: list[str]) -> int:
    status = 0
    for path in paths:
        with open(path, 'rb') as fh:
            data = fh.read()
        fixed = fix(data)
        if fixed != data:
            with open(path, 'wb') as fh:
                fh.write(fixed)
            print(f'fixed end of file: {path}')
            status = 1
    return status


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
