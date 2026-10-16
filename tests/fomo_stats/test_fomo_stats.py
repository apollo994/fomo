"""bin/fomo_stats.py: the gffcompare parser and F1."""

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


def test_parse_accuracy(tmp_path):
    path = tmp_path / 'x.stats'
    path.write_text(STATS)
    acc = fomo_stats.parse_accuracy(str(path))
    assert acc['Base'] == (61.2, 35.4)
    assert acc['Intron_chain'] == (30.0, 15.0)
    assert acc['Locus'] == (32.0, 17.8)


def test_f1():
    assert fomo_stats.f1(0, 0) == 0.0
    assert fomo_stats.f1(50, 50) == 50.0
    assert fomo_stats.f1(30, 15) == pytest.approx(20.0)
