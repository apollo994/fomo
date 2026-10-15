"""bin/filter_gff_by_id.py: drop mode (TD2 coding ids) and keep mode (the filtered-annotation contract)."""

from conftest import gff, run_script, write

ANNOT = gff(
    ('c', 'ens', 'ncRNA_gene', 1, 900, '.', '+', '.', 'ID=gene:G1'),
    ('c', 'ens', 'lnc_RNA', 1, 900, '.', '+', '.', 'ID=transcript:T1;Parent=gene:G1'),
    ('c', 'ens', 'exon', 1, 100, '.', '+', '.', 'Parent=transcript:T1'),
    ('c', 'ens', 'exon', 800, 900, '.', '+', '.', 'Parent=transcript:T1'),
    ('c', 'ens', 'lnc_RNA', 1, 700, '.', '+', '.', 'ID=transcript:T2;Parent=gene:G1'),
    ('c', 'ens', 'exon', 1, 100, '.', '+', '.', 'Parent=transcript:T2'),
    ('c', 'ens', 'exon', 600, 700, '.', '+', '.', 'Parent=transcript:T2'),
    ('c', 'ens', 'ncRNA_gene', 2000, 2900, '.', '-', '.', 'ID=gene:G2'),
    ('c', 'ens', 'lnc_RNA', 2000, 2900, '.', '-', '.', 'ID=transcript:T3;Parent=gene:G2'),
    ('c', 'ens', 'exon', 2000, 2100, '.', '-', '.', 'Parent=transcript:T3'),
    ('c', 'ens', 'exon', 2800, 2900, '.', '-', '.', 'Parent=transcript:T3'),
)


def ids_out(path):
    out = []
    for line in path.read_text().splitlines():
        if line.startswith('#'):
            continue
        attrs = dict(kv.split('=', 1) for kv in line.split('\t')[8].split(';'))
        out.append(attrs.get('ID', '-') + '<' + attrs.get('Parent', ''))
    return out


def test_drop_mode_lineage(tmp_path):
    src = write(tmp_path / 'a.gff3', ANNOT)
    drop = write(tmp_path / 'drop.txt', 'T1\nT3\n')
    out = tmp_path / 'out.gff3'
    r = run_script('filter_gff_by_id.py', '--gff', src, '--drop-ids', drop, '--output', out)
    assert r.returncode == 0, r.stderr
    got = ids_out(out)
    assert 'transcript:T1<gene:G1' not in got and '-<transcript:T1' not in got, 'T1 and its exons dropped'
    assert 'gene:G1<' in got, 'G1 kept: T2 survives'
    assert 'gene:G2<' not in got, 'G2 dropped with its only transcript'
    assert got.count('-<transcript:T2') == 2


def test_drop_mode_empty_list_keeps_everything(tmp_path):
    src = write(tmp_path / 'a.gff3', ANNOT)
    drop = write(tmp_path / 'drop.txt', '')
    out = tmp_path / 'out.gff3'
    run_script('filter_gff_by_id.py', '--gff', src, '--drop-ids', drop, '--output', out, check=True)
    assert out.read_text() == ANNOT


def test_keep_mode_exact(tmp_path):
    src = write(tmp_path / 'a.gff3', ANNOT)
    fasta = write(tmp_path / 'k.fa', '>T2|lnc|S\nACGT\n')
    idmap = write(tmp_path / 'map.tsv', 'T1|lnc|S\ttranscript:T1\nT2|lnc|S\ttranscript:T2\nT3|lnc|S\ttranscript:T3\n')
    out = tmp_path / 'out.gff3'
    r = run_script('filter_gff_by_id.py', '--gff', src, '--keep-fasta', fasta, '--id-map', idmap, '--output', out)
    assert r.returncode == 0, r.stderr
    tx = [x for x in ids_out(out) if x.startswith('transcript:')]
    assert tx == ['transcript:T2<gene:G1'], 'the transcripts of the output are exactly the FASTA records'


def test_keep_mode_fails_on_twice_mapped_header(tmp_path):
    src = write(tmp_path / 'a.gff3', ANNOT)
    fasta = write(tmp_path / 'k.fa', '>T2|lnc|S\nACGT\n')
    idmap = write(tmp_path / 'map.tsv', 'T2|lnc|S\ttranscript:T2\nT2|lnc|S\ttranscript:T1\n')
    r = run_script(
        'filter_gff_by_id.py', '--gff', src, '--keep-fasta', fasta, '--id-map', idmap, '--output', tmp_path / 'o'
    )
    assert r.returncode == 1 and 'maps to both' in r.stderr


def test_keep_mode_fails_on_unmapped_record(tmp_path):
    src = write(tmp_path / 'a.gff3', ANNOT)
    fasta = write(tmp_path / 'k.fa', '>X|lnc|S\nACGT\n')
    idmap = write(tmp_path / 'map.tsv', 'T2|lnc|S\ttranscript:T2\n')
    r = run_script(
        'filter_gff_by_id.py', '--gff', src, '--keep-fasta', fasta, '--id-map', idmap, '--output', tmp_path / 'o'
    )
    assert r.returncode == 1 and 'is not in the id map' in r.stderr


def test_keep_mode_fails_on_id_not_in_gff(tmp_path):
    src = write(tmp_path / 'a.gff3', ANNOT)
    fasta = write(tmp_path / 'k.fa', '>Z|lnc|S\nACGT\n')
    idmap = write(tmp_path / 'map.tsv', 'Z|lnc|S\ttranscript:Z\n')
    r = run_script(
        'filter_gff_by_id.py', '--gff', src, '--keep-fasta', fasta, '--id-map', idmap, '--output', tmp_path / 'o'
    )
    assert r.returncode == 1 and 'have no transcript' in r.stderr


def test_projected_flat_ids(tmp_path):
    rows = gff(
        ('c', 'fomo', 'gene', 1, 900, 60, '+', '.', 'ID=gene-t1|A;source=A'),
        ('c', 'fomo', 'transcript', 1, 900, 60, '+', '.', 'ID=t1|A;Parent=gene-t1|A;source=A'),
        ('c', 'fomo', 'exon', 1, 100, 60, '+', '.', 'ID=t1|A.exon1;Parent=t1|A'),
        ('c', 'fomo', 'gene', 1, 900, 60, '+', '.', 'ID=gene-t1|B;source=B'),
        ('c', 'fomo', 'transcript', 1, 900, 60, '+', '.', 'ID=t1|B;Parent=gene-t1|B;source=B'),
        ('c', 'fomo', 'exon', 1, 100, 60, '+', '.', 'ID=t1|B.exon1;Parent=t1|B'),
    )
    src = write(tmp_path / 'p.gff3', rows)
    drop = write(tmp_path / 'drop.txt', 't1|A\n')
    out = tmp_path / 'out.gff3'
    run_script('filter_gff_by_id.py', '--gff', src, '--drop-ids', drop, '--output', out, check=True)
    text = out.read_text()
    assert 't1|A' not in text and 't1|B' in text, 'source-qualified id drops only that source'


def test_without_drop_list_keeps_everything(tmp_path):
    src = write(tmp_path / 'a.gff3', ANNOT)
    out = tmp_path / 'out.gff3'
    run_script('filter_gff_by_id.py', '--gff', src, '--output', out, check=True)
    assert out.read_text() == ANNOT
