"""bin/fomo_io.py: the shared GFF3 / FASTA / gzip reader, attribute parser, vocabulary and errors."""

import gzip
import subprocess
import sys

import fomo_io
import pytest

from conftest import BIN, write

GFF = (
    '##gff-version 3\n'
    '# a comment\n'
    '\n'
    'c1\tens\tgene\t1\t100\t.\t+\t.\tID=g1\n'
    'not\ta\tfeature\n'
    'c1\tens\tlnc_RNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n'
    '##FASTA\n'
    '>c1\n'
    'ACGT\tA\tB\tC\tD\tE\tF\tG\tH\n'
)


@pytest.mark.parametrize('gz', [False, True])
def test_iter_gff_skips_and_stops_at_fasta(tmp_path, gz):
    path = tmp_path / ('a.gff3.gz' if gz else 'a.gff3')
    if gz:
        with gzip.open(path, 'wt') as fh:
            fh.write(GFF)
    else:
        write(path, GFF)
    rows = list(fomo_io.iter_gff(str(path)))
    assert [n for n, _ in rows] == [4, 6]
    assert rows[1][1][2] == 'lnc_RNA' and rows[1][1][8] == 'ID=t1;Parent=g1'


def test_parse_attrs():
    a = fomo_io.parse_attrs(' ID=t1 ; Parent=g1,g2;Name=a%3Bb;flag;ID=dup; note = with space ;')
    assert a == {'ID': 't1', 'Parent': 'g1,g2', 'Name': 'a%3Bb', 'note': 'with space'}
    assert fomo_io.parse_attrs('.') == {}
    assert fomo_io.parse_attrs('') == {}


def test_split_ids():
    assert fomo_io.split_ids('g1, g2,,') == ['g1', 'g2']
    assert fomo_io.split_ids('') == []


def test_open_text_round_trip(tmp_path):
    for name in ('x.txt', 'x.txt.gz'):
        with fomo_io.open_text(str(tmp_path / name), 'wt') as fh:
            fh.write('a\nb\n')
        with fomo_io.open_text(str(tmp_path / name)) as fh:
            assert fh.read() == 'a\nb\n'
    with open(tmp_path / 'x.txt.gz', 'rb') as fh:
        assert fh.read(2) == b'\x1f\x8b'


@pytest.mark.parametrize('name', ['a.fa', 'a.fa.gz'])
def test_iter_fasta_and_names(tmp_path, name):
    text = 'junk before\n>r1 desc\nAC\nGT\n>r2|lncRNA|sp\n\n>r3\nTT\n'
    path = tmp_path / name
    with fomo_io.open_text(str(path), 'wt') as fh:
        fh.write(text)
    recs = list(fomo_io.iter_fasta(str(path)))
    assert recs == [
        ('r1', '>r1 desc\n', ['AC\n', 'GT\n']),
        ('r2|lncRNA|sp', '>r2|lncRNA|sp\n', ['\n']),
        ('r3', '>r3\n', ['TT\n']),
    ]
    assert fomo_io.fasta_names(str(path)) == ['r1', 'r2|lncRNA|sp', 'r3']


def test_feature_kind():
    assert fomo_io.feature_kind('ncRNA_gene') == 'gene'
    assert fomo_io.feature_kind('NCRNA_GENE') == 'gene'
    assert fomo_io.feature_kind('region') == 'structural'
    assert fomo_io.feature_kind('five_prime_UTR') == 'child'
    assert fomo_io.feature_kind('exon') == 'child'
    for t in ('mRNA', 'lnc_RNA', 'transcript', 'tRNA', 'scRNA', 'misc_RNA'):
        assert fomo_io.is_transcript(t), t


def test_run_maps_fomo_error_to_exit_1(tmp_path):
    script = write(
        tmp_path / 'boom.py',
        f'import sys\nsys.path.insert(0, {str(BIN)!r})\nimport fomo_io\n'
        'def main():\n    fomo_io.die("bad input")\n'
        'fomo_io.run(main)\n',
    )
    r = subprocess.run([sys.executable, str(script)], capture_output=True, text=True)
    assert r.returncode == 1
    assert r.stderr == 'boom: ERROR: bad input\n'


def test_run_passes_the_status(tmp_path):
    script = write(
        tmp_path / 'ok.py',
        f'import sys\nsys.path.insert(0, {str(BIN)!r})\nimport fomo_io\nfomo_io.run(lambda: 3)\n',
    )
    assert subprocess.run([sys.executable, str(script)]).returncode == 3
