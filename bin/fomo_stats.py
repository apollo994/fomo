#!/usr/bin/env python3
"""
Shared helpers for the FOMO bin/ scripts that read gffcompare `.stats` files:
the model-identity filename grammar, the aggregate pseudo-source ids, and the
Sensitivity/Precision parser.

This is a MODULE, not a script — it is imported by `gffcompare_accuracy_mqc.py`
and `run_summary_tables.py`, each of which does

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import fomo_stats

Nextflow only puts `bin/` on PATH, not on Python's import path — but it stages
(or bind-mounts) the whole directory, so a `__file__`-relative insert finds the
siblings on every executor. Do not "simplify" that preamble to a bare import.

The identity regex, AGGREGATE_IDS and the stats parser must agree with each other and
with the pseudo-source id the Nextflow side writes (`allModels_raw`, linked by
GFFCOMPARE_BATCH and GFF_STATS_PROJECTED_BATCH), so they live here once.
"""

import os
import re
from typing import NamedTuple

# Feature levels, in the order gffcompare reports them.
LEVELS = ["Base", "Exon", "Intron", "Intron_chain", "Transcript", "Locus"]

# Pseudo-source ids of aggregate models (every source pooled). They use the per-source
# `from_<id>` filename convention, so per-source tables must skip them.
AGGREGATE_IDS = ("allModels_raw",)

# <target>.from_<source>.<feature_type>[.decoy], after the .gffcompare.stats suffix
# is stripped. Species names contain '_' and digits, hence the non-greedy captures
# anchored on the literal '.from_' and the feature-type alternation.
FNAME_RE = re.compile(r"^(?P<target>.+?)\.from_(?P<source>.+?)\.(?P<ft>lnc_RNA|mRNA)(?P<decoy>\.decoy)?$")


class Identity(NamedTuple):
    """Which model a gffcompare stats file describes."""

    target: str
    source: str  # a real species, or one of AGGREGATE_IDS
    feature_type: str  # 'lnc_RNA' | 'mRNA'
    decoy: bool


def parse_identity(stats_path: str) -> Identity | None:
    """Derive the Identity from a gffcompare stats filename, or None if it does
    not match the grammar."""
    base = os.path.basename(stats_path)
    base = re.sub(r"\.gffcompare\.stats$", "", base)
    base = re.sub(r"\.stats$", "", base)  # tolerate either suffix
    base = re.sub(r"\.gffcompare$", "", base)
    m = FNAME_RE.match(base)
    if not m:
        return None
    return Identity(
        target=m.group("target"),
        source=m.group("source"),
        feature_type=m.group("ft"),
        decoy=bool(m.group("decoy")),
    )


def feature_class(ft: str, decoy: bool) -> str:
    """Presentation label for a (feature_type, decoy) pair."""
    if decoy:
        return "decoy"
    return "lncRNA" if ft == "lnc_RNA" else "mRNA"


def parse_accuracy(stats_path: str) -> dict[str, tuple[float, float]]:
    """Return {level: (sensitivity, precision)} read from a stats file.
    Mirrors MultiQC's native parser: drop '|', normalise 'Intron chain',
    then split — token[0]=level, token[2]=sensitivity, token[3]=precision."""
    out: dict[str, tuple[float, float]] = {}
    with open(stats_path, encoding="utf-8") as fh:
        for line in fh:
            if "level:" not in line:
                continue
            toks = line.replace("|", "").replace("Intron chain", "Intron_chain").split()
            if len(toks) < 4:
                continue
            level = toks[0]
            try:
                sens = float(toks[2])
                prec = float(toks[3])
            except ValueError:
                continue
            out[level] = (sens, prec)
    return out


def transcript_sn_pr(stats_path: str) -> tuple[float, float] | None:
    """(sensitivity, precision) from the Transcript-level row, or None if the file
    has no parseable one. Transcript level is the ranking level throughout FOMO."""
    return parse_accuracy(stats_path).get("Transcript")


def f1(sn: float, pr: float) -> float:
    return 0.0 if (sn + pr) == 0 else 2.0 * sn * pr / (sn + pr)


# Curation report (curate_models.py) and its MultiQC tables (curation_to_mqc.py).
# Locations relative to the reference, in priority order: a gene takes the first that applies.
CURATION_LOCATIONS = (
    'antisense_exonic',
    'intronic',
    'sense_span_overlap',
    'antisense_intronic',
    'divergent',
    'intergenic',
    'no_reference',
)
CURATION_OVERLAPPING = CURATION_LOCATIONS[:4]
# Species buckets of the chain and gene counts: <n>sp for n = --min-species .. this, then ge<n+1>sp.
CURATION_MAX_SPECIES_BUCKET = 4
# Columns of the three MultiQC tables.
CURATION_OBSERVED = (
    'n_models',
    'self_excluded',
    'single_exon',
    'multi_exon_models',
    'species',
    'chains',
    'supported_chains',
    'supported_chains_2sp',
    'supported_chains_3sp',
    'supported_chains_4sp',
    'supported_chains_ge5sp',
    'models_in_supported_chains',
    'genes',
)
CURATION_DISCARDED = (
    'models_removed',
    'models_removed_in_supported_chains',
    'supported_chains_lost',
    'supported_chains_kept_after_removal',
    'genes_lost',
)
CURATION_NOVEL = (
    'genes',
    'transcripts',
    'genes_2sp',
    'genes_3sp',
    'genes_4sp',
    'genes_ge5sp',
    'genes_1chain',
    'genes_2chain',
    'genes_ge3chain',
)
CURATION_LOCATION_COLUMNS = (
    'intergenic',
    'divergent',
    'antisense_exonic',
    'antisense_intronic',
    'intronic',
    'sense_span_overlap',
    'no_reference',
)
CURATION_MEDIANS = ('n_exons', 'transcript_length', 'gene_span', 'intron_length', 'de')
