"""bin/fasta_drop_species.sh: drop one species' records from a `>tid|type|species` FASTA."""

from conftest import run_script, write

FASTA = (
    '>t1|lncRNA|alpha\nACGT\nAC\n'
    '>t2|lncRNA|beta\nGGGG\n'
    '>t3|lncRNA|alpha_2\nTTTT\n'
    '>t4|decoy_lncRNA|alpha\nCCCC\n'
)


def test_drops_only_the_named_species(tmp_path):
    r = run_script('fasta_drop_species.sh', 'alpha', write(tmp_path / 'in.fa', FASTA))
    assert r.returncode == 0, r.stderr
    assert r.stdout == '>t2|lncRNA|beta\nGGGG\n>t3|lncRNA|alpha_2\nTTTT\n'


def test_species_absent_keeps_everything(tmp_path):
    r = run_script('fasta_drop_species.sh', 'gamma', write(tmp_path / 'in.fa', FASTA))
    assert r.returncode == 0, r.stderr
    assert r.stdout == FASTA


def test_every_record_dropped_gives_empty_output(tmp_path):
    r = run_script('fasta_drop_species.sh', 'alpha', write(tmp_path / 'in.fa', '>t1|lncRNA|alpha\nACGT\n'))
    assert r.returncode == 0, r.stderr
    assert r.stdout == ''


def test_usage_error(tmp_path):
    r = run_script('fasta_drop_species.sh', 'alpha')
    assert r.returncode == 1
    assert 'Usage' in r.stderr
