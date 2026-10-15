"""bin/filter_transcript.sh: stranded, spliced and (with --longest) one isoform per gene.

These tests are the behaviour oracle for any reimplementation (plans/27f).
"""

import gzip

import pytest

from conftest import run_script, write

SCRIPT = 'filter_transcript.sh'


def tx(chrom, tid, parent, strand, exons, ftype='lnc_RNA', child_strands=None):
    """Rows of one transcript: the transcript and its exons (optionally with their own strands)."""
    s, e = exons[0][0], exons[-1][1]
    attrs = f'ID={tid}' + (f';Parent={parent}' if parent else '')
    rows = [f'{chrom}\tens\t{ftype}\t{s}\t{e}\t.\t{strand}\t.\t{attrs}']
    for i, (a, b) in enumerate(exons):
        st = child_strands[i] if child_strands else strand
        rows.append(f'{chrom}\tens\texon\t{a}\t{b}\t.\t{st}\t.\tParent={tid}')
    return rows


def gene(chrom, gid, strand, s=1, e=10000, ftype='ncRNA_gene'):
    return [f'{chrom}\tens\t{ftype}\t{s}\t{e}\t.\t{strand}\t.\tID={gid}']


def run(tmp_path, rows, *opts, ttype='lnc_RNA', name='in.gff3'):
    src = write(tmp_path / name, '##gff-version 3\n' + '\n'.join(rows) + '\n')
    counts = tmp_path / 'counts.tsv'
    r = run_script(SCRIPT, *opts, '--counts', counts, '--sample', 'S', src, ttype)
    assert r.returncode == 0, r.stderr
    funnel = counts.read_text().splitlines()
    return r.stdout.splitlines(), funnel


def ids(lines, ftype='lnc_RNA'):
    return [ln.split('\t')[8].split(';')[0][3:] for ln in lines if ln.split('\t')[2] == ftype]


def test_spliced_filter_and_funnel(tmp_path):
    rows = (
        gene('c', 'G1', '+') + tx('c', 'T1', 'G1', '+', [(1, 100), (200, 300)]) + tx('c', 'T2', 'G1', '+', [(1, 300)])
    )
    out, funnel = run(tmp_path, rows)
    assert ids(out) == ['T1'], 'mono-exonic T2 dropped'
    assert funnel == ['Sample\tfeature\tn_transcripts\tn_stranded\tn_spliced\tn_selected', 'S\tlnc_RNA\t2\t2\t1\t1']
    assert out[0].split('\t')[2] == 'ncRNA_gene', 'parent gene printed'
    assert not any(ln.startswith('#') for ln in out), 'no header or comment lines are emitted'


def test_unstranded_and_mixed_strand_dropped(tmp_path):
    rows = (
        gene('c', 'G1', '+')
        + tx('c', 'Tq', 'G1', '?', [(1, 100), (200, 300)])
        + tx('c', 'Tmix', 'G1', '+', [(1, 100), (200, 300)], child_strands=['+', '-'])
        + tx('c', 'Tdot', 'G1', '+', [(1, 100), (200, 300)], child_strands=['+', '.'])
        + gene('c', 'Gbad', '?')
        + tx('c', 'Tbadparent', 'Gbad', '+', [(1, 100), (200, 300)])
        + tx('c', 'Tok', 'G1', '+', [(1, 100), (200, 300)])
    )
    out, funnel = run(tmp_path, rows)
    assert ids(out) == ['Tok']
    assert funnel[1] == 'S\tlnc_RNA\t5\t1\t1\t1'


def test_longest_isoform_by_spliced_length(tmp_path):
    rows = (
        gene('c', 'G1', '+')
        + tx('c', 'Tshort_span_long_exons', 'G1', '+', [(1, 400), (500, 900)])
        + tx('c', 'Tlong_span', 'G1', '+', [(1, 100), (8000, 8100)])
    )
    out, funnel = run(tmp_path, rows, '--longest')
    assert ids(out) == ['Tshort_span_long_exons'], 'summed exon length decides, not the genomic span'
    assert funnel[1] == 'S\tlnc_RNA\t2\t2\t2\t1'


def test_longest_tie_goes_to_smallest_id_bytewise(tmp_path):
    rows = (
        gene('c', 'G1', '+')
        + tx('c', 'b', 'G1', '+', [(1, 100), (200, 300)])
        + tx('c', 'B', 'G1', '+', [(1, 100), (200, 300)])
        + tx('c', 'a', 'G1', '+', [(1, 100), (200, 300)])
    )
    out, _ = run(tmp_path, rows, '--longest')
    assert ids(out) == ['B'], "byte order: 'B' < 'a' < 'b'"


def test_longest_ignores_monoexonic_longest(tmp_path):
    rows = (
        gene('c', 'G1', '+')
        + tx('c', 'Tmono', 'G1', '+', [(1, 5000)])
        + tx('c', 'Tspliced', 'G1', '+', [(1, 100), (200, 300)])
    )
    out, _ = run(tmp_path, rows, '--longest')
    assert ids(out) == ['Tspliced'], 'the longest spliced isoform survives a longer mono-exonic one'


def test_one_per_gene_and_parentless_transcripts(tmp_path):
    rows = (
        gene('c', 'G1', '+')
        + gene('c', 'G2', '-')
        + tx('c', 'A1', 'G1', '+', [(1, 100), (200, 300)])
        + tx('c', 'A2', 'G1', '+', [(1, 100), (200, 400)])
        + tx('c', 'B1', 'G2', '-', [(1, 100), (200, 300)])
        + tx('c', 'Orphan1', '', '+', [(1, 100), (200, 300)])
        + tx('c', 'Orphan2', '', '+', [(1, 100), (200, 300)])
    )
    out, funnel = run(tmp_path, rows, '--longest')
    assert sorted(ids(out)) == ['A2', 'B1', 'Orphan1', 'Orphan2'], (
        'a transcript without Parent competes only with itself'
    )
    assert funnel[1] == 'S\tlnc_RNA\t5\t5\t5\t4'


def test_shared_exon_credits_every_parent(tmp_path):
    rows = gene('c', 'G1', '+') + [
        'c\tens\tlnc_RNA\t1\t300\t.\t+\t.\tID=T1;Parent=G1',
        'c\tens\tlnc_RNA\t1\t300\t.\t+\t.\tID=T2;Parent=G1',
        'c\tens\texon\t1\t100\t.\t+\t.\tParent=T1,T2',
        'c\tens\texon\t200\t300\t.\t+\t.\tParent=T1,T2',
    ]
    out, funnel = run(tmp_path, rows)
    assert ids(out) == ['T1', 'T2'] and funnel[1] == 'S\tlnc_RNA\t2\t2\t2\t2'
    assert sum(1 for ln in out if '\texon\t' in ln) == 2, 'shared exon rows printed once'


def test_without_longest_keeps_every_spliced_isoform(tmp_path):
    rows = (
        gene('c', 'G1', '+')
        + tx('c', 'A1', 'G1', '+', [(1, 100), (200, 300)])
        + tx('c', 'A2', 'G1', '+', [(1, 100), (200, 400)])
    )
    out, funnel = run(tmp_path, rows)
    assert ids(out) == ['A1', 'A2'] and funnel[1] == 'S\tlnc_RNA\t2\t2\t2\t2'


def test_other_types_and_cds_children(tmp_path):
    rows = (
        gene('c', 'G1', '+', ftype='gene')
        + tx('c', 'M1', 'G1', '+', [(1, 100), (200, 300)], ftype='mRNA')
        + ['c\tens\tCDS\t50\t100\t.\t+\t0\tParent=M1']
        + tx('c', 'L1', 'G1', '+', [(1, 100), (200, 300)])
    )
    out, _ = run(tmp_path, rows, ttype='mRNA')
    assert ids(out, 'mRNA') == ['M1'] and not ids(out, 'lnc_RNA')
    assert any('\tCDS\t' in ln for ln in out), 'CDS children of a selected transcript are kept'


def test_empty_input_reports_zero(tmp_path):
    src = write(tmp_path / 'e.gff3', '')
    counts = tmp_path / 'c.tsv'
    r = run_script(SCRIPT, '--counts', counts, '--sample', 'S', src, 'lnc_RNA')
    assert r.returncode == 0 and r.stdout == ''
    assert counts.read_text().splitlines()[1] == 'S\tlnc_RNA\t0\t0\t0\t0'


def test_output_order_is_input_order(tmp_path):
    rows = (
        gene('c', 'G1', '+')
        + tx('c', 'T1', 'G1', '+', [(1, 100), (200, 300)])
        + gene('c', 'G2', '+')
        + tx('c', 'T2', 'G2', '+', [(1, 100), (200, 300)])
    )
    out, _ = run(tmp_path, rows)
    assert out == rows


@pytest.mark.parametrize('opts', [[], ['--longest']])
def test_real_annotation_slice_is_stable(tmp_path, opts):
    """A slice of a real Ensembl annotation (committed test data) runs and selects a subset."""
    from conftest import BIN

    src = BIN.parent / 'assets/test_data/Lycaena_alciphron_282377/GCA_964273875.1/Lycaena_alciphron_282377.gff3.gz'
    plain = tmp_path / 'a.gff3'
    with gzip.open(src, 'rt') as fh:
        plain.write_text(fh.read())
    out, funnel = run(tmp_path, plain.read_text().splitlines(), *opts, name='real.gff3')
    n = [int(x) for x in funnel[1].split('\t')[2:]]
    assert n[0] >= n[1] >= n[2] >= n[3] > 0
    assert len(ids(out)) == n[3]
