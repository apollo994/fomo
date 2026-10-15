"""fomo_stats record helpers (fomo-data/1) and the gffcompare stats parser."""

import json
from pathlib import Path

import fomo_io
import fomo_stats
import pytest

FIXTURES = Path(__file__).parent / 'fixtures'


def test_round_trip(tmp_path):
    a = fomo_stats.make_record('curation', {'target': 'T', 'side': 'target', 'track': 'lnc_RNA'}, {'x': 1})
    b = fomo_stats.make_record('fasta_stats', {'species': 'S', 'side': 'source', 'track': 'mRNA'}, {'n': 2})
    fomo_stats.write_records(tmp_path / 'one.data.json', a)
    (tmp_path / 'sub').mkdir()
    fomo_stats.write_records(tmp_path / 'sub' / 'two.data.json', [b])
    (tmp_path / 'ignored.json').write_text('{}')
    assert fomo_stats.read_records(tmp_path) == [a, b]


def test_missing_identity_key_fails():
    with pytest.raises(fomo_io.FomoError, match='missing identity'):
        fomo_stats.make_record('accuracy', {'target': 'T', 'side': 'target', 'track': 'lnc_RNA'}, {})


def test_unknown_section_and_key_fail(tmp_path):
    with pytest.raises(fomo_io.FomoError, match='unknown record section'):
        fomo_stats.make_record('nope', {}, {})
    with pytest.raises(fomo_io.FomoError, match='unknown identity key'):
        fomo_stats.make_record('curation', {'target': 'T', 'side': 'target', 'track': 't', 'colour': 'x'}, {})
    (tmp_path / 'bad.data.json').write_text(json.dumps({'schema': 'other/1'}))
    with pytest.raises(fomo_io.FomoError, match='not a fomo-data/1 record'):
        fomo_stats.read_records(tmp_path)


def test_parse_gffcompare_stats():
    acc = fomo_stats.parse_gffcompare_stats(
        str(FIXTURES / 'Lycaena_alciphron_282377.from_allModels_raw.lnc_RNA.gffcompare.stats')
    )
    assert acc['levels']['Transcript'] == {'sn': 28.8, 'pr': 28.3}
    assert list(acc['levels']) == fomo_stats.LEVELS
    assert acc['counts']['query_mrnas'] == 60 and acc['counts']['reference_loci'] == 59
    assert acc['counts']['matching_intron_chains'] == 17
    assert (acc['counts']['missed_exons'], acc['counts']['missed_exons_total']) == (84, 131)


def test_identity_arg():
    assert fomo_stats.identity_arg('{"species": "S"}') == {'species': 'S'}
    with pytest.raises(fomo_io.FomoError):
        fomo_stats.identity_arg('[1]')
