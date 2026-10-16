"""Shared helpers for the bin/ script tests."""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parent.parent / 'bin'
sys.path.insert(0, str(BIN))


def run_script(name, *args, input=None, check=False):
    """Run bin/<name> (Python scripts with this interpreter) and return the CompletedProcess."""
    path = BIN / name
    cmd = [sys.executable, str(path)] if name.endswith('.py') else [str(path)]
    return subprocess.run(cmd + [str(a) for a in args], capture_output=True, text=True, input=input, check=check)


def needs(*tools):
    """Skip marker for tests that need external tools on PATH."""
    missing = [t for t in tools if shutil.which(t) is None]
    return pytest.mark.skipif(bool(missing), reason=f"needs {', '.join(missing)} on PATH")


def write(path, text):
    path = Path(path)
    path.write_text(text, encoding='utf-8')
    return path


def gff(*rows):
    """GFF3 text from rows given as 9-tuples (or tab-joined strings)."""
    lines = ['##gff-version 3']
    lines += ['\t'.join(map(str, r)) if not isinstance(r, str) else r for r in rows]
    return '\n'.join(lines) + '\n'
