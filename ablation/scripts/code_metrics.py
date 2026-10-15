#!/usr/bin/env python3
"""Size, comment density and "documentation debt" markers for every tracked source file of the
pipeline (ablation/CODE_REVIEW.md). Run from the repository root.

Columns
  code / comment / blank   line counts (a line holding only a comment is a comment line; a
                           Python docstring counts as comment)
  ratio                    comment lines per code line
  plans                    references to design plans ("plans/NN")
  history                  phrases that narrate past states rather than current behaviour
  incident                 references to specific runs, jobs, clades or measurements

Usage: code_metrics.py > code_metrics.tsv
"""
import ast
import re
import subprocess
import sys

HISTORY = re.compile(r'\b(used to|no longer|was (?:removed|replaced|dropped)|since plans|pre-plans|previously|'
                     r'the old |old [A-Z_]{4,}|before plans|went with|is gone|are gone|anymore)\b', re.I)
INCIDENT = re.compile(r'(seen for real|hit for real|verified|jobid|job \d{6,}|\b\d{8}\b|agitated_|thirsty_|'
                      r'lepidoptera|lycaenidae|mollusca|apiaceae|acari|drosophil|cephalopod|vampyroteuthis|'
                      r'daucus|tick chromosome|measured|observed:|trace[sd]? of|\brun `)', re.I)
PLANS = re.compile(r'plans/\d+')

LANG = {'.nf': '//', '.config': '//', '.groovy': '//', '.sh': '#', '.py': '#', '.yml': '#', '.yaml': '#'}


def docstring_lines(src):
    """Line numbers covered by Python docstrings (module, class, function)."""
    out = set()
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return out
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(getattr(first, 'value', None), ast.Constant) \
                    and isinstance(first.value.value, str):
                out.update(range(first.lineno, first.end_lineno + 1))
    return out


def classify(path):
    ext = '.' + path.rsplit('.', 1)[-1] if '.' in path else ''
    if path == 'main.nf' or path.endswith('nextflow.config'):
        ext = '.nf'
    mark = LANG.get(ext)
    if mark is None:
        return None
    src = open(path, encoding='utf-8', errors='replace').read()
    lines = src.split('\n')
    doc = docstring_lines(src) if ext == '.py' else set()
    code = comment = blank = 0
    in_block = False
    comment_text = []
    for i, line in enumerate(lines, 1):
        s = line.strip()
        if not s:
            blank += 1
            continue
        if i in doc:
            comment += 1
            comment_text.append(s)
            continue
        if mark == '//':
            if in_block:
                comment += 1
                comment_text.append(s)
                in_block = '*/' not in s
                continue
            if s.startswith('/*'):
                comment += 1
                comment_text.append(s)
                in_block = '*/' not in s
                continue
        if s.startswith(mark) and not s.startswith('#!'):
            comment += 1
            comment_text.append(s)
            continue
        code += 1
        # trailing comment
        idx = line.find(' ' + mark + ' ')
        if idx > 0:
            comment_text.append(line[idx:])
    text = '\n'.join(comment_text)
    return dict(code=code, comment=comment, blank=blank,
                ratio=comment / code if code else float('inf'),
                plans=len(PLANS.findall(text)), history=len(HISTORY.findall(text)),
                incident=len(INCIDENT.findall(text)))


def main():
    files = subprocess.run(['git', 'ls-files'], capture_output=True, text=True, check=True).stdout.split()
    files = [f for f in files if not f.startswith(('plans/', 'assets/test_data/', 'ablation/'))
             and not f.endswith('gff-feature-stats')]
    rows = []
    for f in files:
        m = classify(f)
        if m:
            rows.append((f, m))
    print('file\tcode\tcomment\tblank\tratio\tplans\thistory\tincident')
    tot = dict(code=0, comment=0, plans=0, history=0, incident=0)
    for f, m in sorted(rows, key=lambda x: -x[1]['ratio']):
        print(f"{f}\t{m['code']}\t{m['comment']}\t{m['blank']}\t{m['ratio']:.2f}\t{m['plans']}\t{m['history']}\t{m['incident']}")
        for k in tot:
            tot[k] += m[k]
    print(f"TOTAL\t{tot['code']}\t{tot['comment']}\t\t{tot['comment'] / tot['code']:.2f}\t{tot['plans']}\t"
          f"{tot['history']}\t{tot['incident']}", file=sys.stderr)


if __name__ == '__main__':
    main()
