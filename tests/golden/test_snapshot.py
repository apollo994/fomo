"""tests/golden/snapshot.py: the normalisation rules of the golden gate."""

import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location('snapshot', Path(__file__).with_name('snapshot.py'))
snapshot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(snapshot)

LOCI = 'XLOC_000201\tc1[+]10-90\tgene:g1|transcript:t1\tq2|sp,q1|sp\nXLOC_000202\tc1[-]100-200\t-\tq3|sp\n'


def write(tmp_path, name, text):
    path = tmp_path / 'gffcompare' / name
    path.parent.mkdir(exist_ok=True)
    path.write_text(text)
    return path


def test_loci_tie_order_is_ignored(tmp_path):
    a = write(tmp_path, 'a.loci', LOCI)
    b = write(tmp_path, 'b.loci', LOCI.replace('q2|sp,q1|sp', 'q1|sp,q2|sp'))
    assert snapshot.normalise(a) == snapshot.normalise(b)
    assert snapshot.normalise(a)[0] == 'sorted_entries'


def test_loci_content_still_counts(tmp_path):
    a = write(tmp_path, 'a.loci', LOCI)
    for changed in (LOCI.replace('q2|sp,', ''), LOCI.replace('XLOC_000202', 'XLOC_000203'), LOCI.replace('t1', 't9')):
        assert snapshot.normalise(write(tmp_path, 'b.loci', changed)) != snapshot.normalise(a)


def test_other_gffcompare_files_stay_bytes(tmp_path):
    tracking = write(tmp_path, 'a.tracking', 'TCONS_1\tXLOC_1\tg1|t1\t=\tq:q1|sp\n')
    assert snapshot.normalise(tracking)[0] == 'bytes'
