"""bin/fomo_stats.py: the stats filename grammar, the gffcompare parser and F1."""

import pytest

import fomo_stats

STATS = """# gffcompare v0.12.6 | Command line was:
#gffcompare -M --no-merge -r ref.gff3 -o T.from_S.lnc_RNA.gffcompare q.gff3
#

#= Summary for dataset: q.gff3
#     Query mRNAs :     100 in      90 loci  (100 multi-exon transcripts)
#            (  5 multi-transcript loci, ~1.1 transcripts per locus)
# Reference mRNAs :      50 in      50 loci  (50 multi-exon)
# Super-loci w/ reference transcripts:       40
#-----------------| Sensitivity | Precision  |
        Base level:    61.2     |    35.4    |
        Exon level:    45.0     |    30.1    |
      Intron level:    52.3     |    40.2    |
Intron chain level:    30.0     |    15.0    |
  Transcript level:    30.0     |    15.0    |
       Locus level:    32.0     |    17.8    |
"""


@pytest.mark.parametrize(
    'name, expected',
    [
        ('T_1.from_S_2.lnc_RNA.gffcompare.stats', ('T_1', 'S_2', 'lnc_RNA', False)),
        ('T_1.from_S_2.mRNA.decoy.gffcompare.stats', ('T_1', 'S_2', 'mRNA', True)),
        ('T_1.from_allModels_raw.lnc_RNA.gffcompare.stats', ('T_1', 'allModels_raw', 'lnc_RNA', False)),
        ('T_1.from_S_2.lnc_RNA.stats', ('T_1', 'S_2', 'lnc_RNA', False)),
    ],
)
def test_parse_identity(name, expected):
    ident = fomo_stats.parse_identity(f'/some/dir/{name}')
    assert (ident.target, ident.source, ident.feature_type, ident.decoy) == expected


def test_parse_identity_aggregate_and_self():
    assert fomo_stats.parse_identity('T.from_allModels_raw.lnc_RNA.stats').is_aggregate
    assert fomo_stats.parse_identity('T.from_T.lnc_RNA.stats').is_self


@pytest.mark.parametrize('name', ['T.lnc_RNA.stats', 'T.from_S.snoRNA.stats', 'garbage'])
def test_parse_identity_rejects_other_names(name):
    assert fomo_stats.parse_identity(name) is None


def test_parse_accuracy(tmp_path):
    path = tmp_path / 'x.stats'
    path.write_text(STATS)
    acc = fomo_stats.parse_accuracy(str(path))
    assert acc['Base'] == (61.2, 35.4)
    assert acc['Intron_chain'] == (30.0, 15.0)
    assert acc['Locus'] == (32.0, 17.8)
    assert fomo_stats.transcript_sn_pr(str(path)) == (30.0, 15.0)


def test_f1():
    assert fomo_stats.f1(0, 0) == 0.0
    assert fomo_stats.f1(50, 50) == 50.0
    assert fomo_stats.f1(30, 15) == pytest.approx(20.0)


@pytest.mark.parametrize(
    'ft, decoy, label',
    [('lnc_RNA', False, 'lncRNA'), ('mRNA', False, 'mRNA'), ('lnc_RNA', True, 'decoy'), ('mRNA', True, 'decoy')],
)
def test_feature_class(ft, decoy, label):
    assert fomo_stats.feature_class(ft, decoy) == label
