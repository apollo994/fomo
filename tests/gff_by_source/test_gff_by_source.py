"""bin/gff_by_source.py: split a merged projection GFF3 by its source= attribute."""

from conftest import gff, run_script, write

TEMPLATE = 'T.from_{source}.lnc_RNA.projected.gff3'


def model(tid, src, chrom='chr1', start=100):
    uid = f'{tid}|{src}'
    return [
        (chrom, 'fomo', 'gene', start, start + 300, 60, '+', '.', f'ID=gene-{uid};source={src};gene_class=lnc_RNA'),
        (
            chrom,
            'fomo',
            'transcript',
            start,
            start + 300,
            60,
            '+',
            '.',
            f'ID={uid};Parent=gene-{uid};source={src};gene_class=lnc_RNA',
        ),
        (chrom, 'fomo', 'exon', start, start + 100, 60, '+', '.', f'ID={uid}.exon1;Parent={uid}'),
        (chrom, 'fomo', 'exon', start + 200, start + 300, 60, '+', '.', f'ID={uid}.exon2;Parent={uid}'),
    ]


def test_split_per_source(tmp_path):
    rows = model('t1', 'A') + model('t2', 'B', start=1000) + model('t3', 'A', start=2000)
    src = write(tmp_path / 'in.gff3', gff(*rows))
    out = tmp_path / 'out'
    r = run_script('gff_by_source.py', '--gff', src, '--outdir', out, '--name-template', TEMPLATE)
    assert r.returncode == 0, r.stderr
    assert sorted(p.name for p in out.iterdir()) == [
        'T.from_A.lnc_RNA.projected.gff3',
        'T.from_B.lnc_RNA.projected.gff3',
    ]
    a = (out / 'T.from_A.lnc_RNA.projected.gff3').read_text().splitlines()
    assert a[0] == '##gff-version 3'
    assert len(a) == 1 + 8, 'two models of A: gene, transcript, two exons each'
    assert all('|A' in line for line in a[1:]), 'children inherit the source through Parent'
    b = (out / 'T.from_B.lnc_RNA.projected.gff3').read_text().splitlines()
    assert [line.split('\t')[2] for line in b[1:]] == ['gene', 'transcript', 'exon', 'exon']


def test_records_keep_input_order(tmp_path):
    rows = model('t1', 'A', start=100) + model('t2', 'A', start=500)
    src = write(tmp_path / 'in.gff3', gff(*rows))
    out = tmp_path / 'out'
    run_script('gff_by_source.py', '--gff', src, '--outdir', out, '--name-template', TEMPLATE, check=True)
    lines = (out / 'T.from_A.lnc_RNA.projected.gff3').read_text().splitlines()[1:]
    assert lines == ['\t'.join(map(str, r)) for r in rows]


def test_duplicate_transcript_id_fails(tmp_path):
    rows = model('t1', 'A')
    dup = [r for r in model('t1', 'B', start=900)]
    dup[1] = dup[1][:8] + ('ID=t1|A;Parent=gene-t1|B;source=B;gene_class=lnc_RNA',)
    src = write(tmp_path / 'in.gff3', gff(*(rows + dup)))
    r = run_script('gff_by_source.py', '--gff', src, '--outdir', tmp_path / 'out', '--name-template', TEMPLATE)
    assert r.returncode == 1
    assert "duplicate transcript ID 't1|A'" in r.stderr and "'A'" in r.stderr and "'B'" in r.stderr


def test_duplicate_transcript_id_from_one_source_fails(tmp_path):
    rows = model('t1', 'A') + model('t1', 'A', start=900)
    src = write(tmp_path / 'in.gff3', gff(*rows))
    r = run_script('gff_by_source.py', '--gff', src, '--outdir', tmp_path / 'out', '--name-template', TEMPLATE)
    assert r.returncode == 1 and "duplicate transcript ID 't1|A'" in r.stderr


def test_empty_input_writes_nothing(tmp_path):
    src = write(tmp_path / 'in.gff3', '##gff-version 3\n')
    out = tmp_path / 'out'
    r = run_script('gff_by_source.py', '--gff', src, '--outdir', out, '--name-template', TEMPLATE)
    assert r.returncode == 0
    assert list(out.iterdir()) == []


def test_unattributable_record_fails(tmp_path):
    src = write(tmp_path / 'in.gff3', gff(('chr1', 'fomo', 'exon', 1, 10, '.', '+', '.', 'ID=e;Parent=unknown')))
    r = run_script('gff_by_source.py', '--gff', src, '--outdir', tmp_path / 'out', '--name-template', TEMPLATE)
    assert r.returncode == 1 and 'cannot attribute' in r.stderr


def test_template_needs_placeholder(tmp_path):
    src = write(tmp_path / 'in.gff3', gff(*model('t1', 'A')))
    r = run_script('gff_by_source.py', '--gff', src, '--outdir', tmp_path / 'out', '--name-template', 'x.gff3')
    assert r.returncode != 0


def test_manifest_lists_every_output_and_the_aggregate(tmp_path):
    rows = model('t1', 'B') + model('t2', 'A', start=900)
    src = write(tmp_path / 'T.allModels.lnc_RNA.raw.gff3', gff(*rows))
    out, manifest = tmp_path / 'out', tmp_path / 'm.tsv'
    args = ['--gff', src, '--outdir', out, '--name-template', TEMPLATE, '--manifest', manifest]
    run_script('gff_by_source.py', *args, '--aggregate-id', 'allModels_raw', check=True)
    lines = [ln.split('\t') for ln in manifest.read_text().splitlines()]
    assert lines[0] == ['file', 'source']
    assert lines[1] == ['T.allModels.lnc_RNA.raw.gff3', 'allModels_raw']
    assert lines[2:] == [[TEMPLATE.replace('{source}', s), s] for s in ('A', 'B')]
    assert all((out / f).exists() for f, _ in lines[2:])


def test_manifest_of_an_empty_input_is_header_only(tmp_path):
    src = write(tmp_path / 'in.gff3', '##gff-version 3\n')
    manifest = tmp_path / 'm.tsv'
    run_script(
        'gff_by_source.py',
        '--gff',
        src,
        '--outdir',
        tmp_path / 'o',
        '--name-template',
        TEMPLATE,
        '--manifest',
        manifest,
        check=True,
    )
    assert manifest.read_text() == 'file\tsource\n'
