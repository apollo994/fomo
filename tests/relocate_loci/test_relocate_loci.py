"""bin/relocate_loci.py: decoys keep a model's shape and move into intergenic space, reproducibly."""

from conftest import run_script, write

MODELS = (
    '\n'.join(
        row
        for i, (start, exons) in enumerate(
            [
                (100, [(100, 200), (400, 500)]),
                (1000, [(1000, 1050), (1200, 1300), (1500, 1550)]),
                (3000, [(3000, 3100)]),
            ]
        )
        for row in [f'c\tens\tlnc_RNA\t{start}\t{exons[-1][1]}\t.\t+\t.\tID=t{i};Parent=g{i}']
        + [f'c\tens\texon\t{a}\t{b}\t.\t+\t.\tParent=t{i}' for a, b in exons]
    )
    + '\n'
)
BED = 'c\t10000\t20000\nc\t30000\t31000\nd\t0\t50000\n'


def run(tmp_path, *extra, seed=7, name='out.gff3'):
    gff, bed = write(tmp_path / 'm.gff3', '##gff-version 3\n' + MODELS), write(tmp_path / 'i.bed', BED)
    out = tmp_path / name
    r = run_script(
        'relocate_loci.py', '--input-gff', gff, '--intergenic-bed', bed, '--output-gff', out, '--seed', seed, *extra
    )
    assert r.returncode == 0, r.stderr
    return out.read_text()


def models(text):
    tx, exons = {}, {}
    for line in text.splitlines():
        if line.startswith('#'):
            continue
        f = line.split('\t')
        attrs = dict(kv.split('=', 1) for kv in f[8].split(';'))
        if f[2] == 'lnc_RNA':
            tx[attrs['ID']] = f
        elif f[2] == 'exon':
            exons.setdefault(attrs['Parent'], []).append((int(f[3]), int(f[4]), f[0]))
    return tx, exons


def test_same_seed_same_output(tmp_path):
    assert run(tmp_path, name='a.gff3') == run(tmp_path, name='b.gff3')


def test_shape_kept_and_inside_intergenic(tmp_path):
    tx, exons = models(run(tmp_path))
    assert set(tx) == {'t0.reloc', 't1.reloc'}, 'single-exon models are not relocated; ids get .reloc'
    shapes = {'t0.reloc': [(0, 100), (300, 400)], 't1.reloc': [(0, 50), (200, 300), (500, 550)]}
    intervals = [('c', 10001, 20000), ('c', 30001, 31000), ('d', 1, 50000)]
    for tid, ex in exons.items():
        s0 = ex[0][0]
        assert [(a - s0, b - s0) for a, b, _ in ex] == shapes[tid], 'exon and intron lengths preserved'
        chrom = ex[0][2]
        assert any(c == chrom and lo <= ex[0][0] and ex[-1][1] <= hi for c, lo, hi in intervals), 'inside one interval'


def test_output_grouped_by_seqid(tmp_path):
    seqids = [ln.split('\t')[0] for ln in run(tmp_path).splitlines() if not ln.startswith('#')]
    blocks = [s for i, s in enumerate(seqids) if i == 0 or s != seqids[i - 1]]
    assert len(blocks) == len(set(blocks))


def test_cap(tmp_path):
    tx, _ = models(run(tmp_path, '--decoy-cap', 1))
    assert len(tx) == 1
