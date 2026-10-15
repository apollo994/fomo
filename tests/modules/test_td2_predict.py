"""TD2_PREDICT treats "nothing to filter" as a result, not a failure: an empty FASTA and a FASTA
whose transcripts have no complete ORF both give an empty .pep without running TD2.Predict (which
dies when PSAURON has nothing to score)."""

import re
import subprocess
from pathlib import Path

MODULE = Path(__file__).resolve().parents[2] / 'modules/local/td2_predict.nf'
LONGORFS_WITH_ORFS = 'TD2.LongOrfs() { mkdir -p td2_work; printf ">t1.p1\\nMK*\\n" > td2_work/longest_orfs.pep; }\n'
LONGORFS_WITHOUT = 'TD2.LongOrfs() { mkdir -p td2_work; : > td2_work/longest_orfs.pep; }\n'
PREDICT = 'TD2.Predict() { echo called > predict_called; printf ">t1.p1\\nMK*\\n" > in.fa.TD2.pep; }\n'


def script():
    body = re.search(r'\n    script:\n.*?"""\n(.*?)\n    """', MODULE.read_text(), re.S).group(1)
    body = body.replace('${fasta}', 'in.fa').replace('${longorfs_args}', '').replace('${predict_args}', '')
    assert not re.search(r'(?<!\\)\$\{[a-z]', body)
    return body.replace('\\$', '$')


def run(tmp_path, fasta, stubs):
    (tmp_path / 'in.fa').write_text(fasta)
    r = subprocess.run(['bash', '-c', stubs + script()], cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return (tmp_path / 'in.fa.TD2.pep').read_text(), (tmp_path / 'predict_called').exists()


def test_empty_fasta_skips_td2(tmp_path):
    assert run(tmp_path, '', 'TD2.LongOrfs() { exit 9; }\n' + PREDICT) == ('', False)


def test_no_complete_orf_skips_predict(tmp_path):
    assert run(tmp_path, '>t1\nACGT\n', LONGORFS_WITHOUT + PREDICT) == ('', False)


def test_orfs_run_predict(tmp_path):
    pep, called = run(tmp_path, '>t1\nACGT\n', LONGORFS_WITH_ORFS + PREDICT)
    assert called and pep.startswith('>t1.p1')
