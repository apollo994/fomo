#!/usr/bin/env python3
"""
Shared statistics helpers of the FOMO bin/ scripts: the data-record format (fomo-data/1), the
gffcompare `.stats` parser, the model-identity filename grammar, and the curation constants.

A module, not a script. Importers put bin/ on the import path first:

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import fomo_stats

Nextflow puts bin/ on PATH but not on Python's import path; it stages or mounts the whole
directory, so a `__file__`-relative insert finds the siblings on every executor.
"""

import json
import os
import re
import sys
from pathlib import Path
from typing import NamedTuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fomo_io  # noqa: E402

# ── data records ──────────────────────────────────────────────────────────────────
# A record is {"schema", "section", "identity", "data"}: what was measured (identity keys drawn
# from IDENTITY_KEYS) and the measurement. A *.data.json file holds one record or a list.
SCHEMA = 'fomo-data/1'
IDENTITY_KEYS = ('species', 'side', 'target', 'source', 'track', 'label', 'order', 'kind', 'stage')
# section -> identity keys every record of it carries
SECTIONS = {
    'transcript_filter': ('species', 'side', 'track'),
    'gffstats': ('species', 'side', 'kind'),
    'td2_coding': ('species', 'side', 'track', 'stage'),
    'alignment': ('target', 'side', 'track'),
    'fasta_stats': ('species', 'side', 'track'),
    'accuracy': ('target', 'source', 'side', 'track'),
    'curation': ('target', 'side', 'track'),
}


def make_record(section: str, identity: dict, data) -> dict:
    """A validated record: a known section, its required identity keys, no unknown key."""
    if section not in SECTIONS:
        fomo_io.die(f'unknown record section {section!r}')
    unknown = sorted(set(identity) - set(IDENTITY_KEYS))
    if unknown:
        fomo_io.die(f'{section} record: unknown identity key(s) {unknown}')
    missing = [k for k in SECTIONS[section] if identity.get(k) in (None, '')]
    if missing:
        fomo_io.die(f'{section} record: missing identity key(s) {missing} in {identity}')
    return {'schema': SCHEMA, 'section': section, 'identity': dict(identity), 'data': data}


def write_records(path, records) -> None:
    """Write one record, or a list of records, as JSON."""
    with fomo_io.open_text(str(path), 'wt') as out:
        json.dump(records, out, indent=1, sort_keys=False)
        out.write('\n')


def read_records(directory) -> list[dict]:
    """Every record of every *.data.json under directory, each validated."""
    records = []
    for path in sorted(Path(directory).rglob('*.data.json')):
        with fomo_io.open_text(str(path)) as fh:
            doc = json.load(fh)
        for rec in doc if isinstance(doc, list) else [doc]:
            if not isinstance(rec, dict) or rec.get('schema') != SCHEMA:
                fomo_io.die(f'{path}: not a {SCHEMA} record')
            records.append(make_record(rec.get('section'), rec.get('identity') or {}, rec.get('data')))
    return records


def identity_arg(text: str) -> dict:
    """The --identity argument of the record-writing scripts: a JSON object."""
    try:
        identity = json.loads(text)
    except json.JSONDecodeError as exc:
        fomo_io.die(f'--identity is not JSON: {exc}')
    if not isinstance(identity, dict):
        fomo_io.die('--identity must be a JSON object')
    return identity


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


GFFCOMPARE_COUNT_RE = re.compile(r'^#?\s*(Query mRNAs|Reference mRNAs)\s*:\s*(\d+)\s+in\s+(\d+)\s+loci')
GFFCOMPARE_MATCH_RE = re.compile(r'^\s*Matching (intron chains|transcripts|loci):\s*(\d+)')
GFFCOMPARE_MISSED_RE = re.compile(r'^\s*(Missed|Novel) (exons|introns|loci):\s*(\d+)/(\d+)')


def parse_gffcompare_stats(stats_path: str) -> dict:
    """The accuracy payload of a gffcompare `.stats` file: Sn/Pr per level, the query and
    reference counts, the matching counts, and the missed/novel counts that are present."""
    counts: dict = {}
    with open(stats_path, encoding='utf-8') as fh:
        for line in fh:
            m = GFFCOMPARE_COUNT_RE.match(line)
            if m:
                key = 'query' if m.group(1).startswith('Query') else 'reference'
                counts[f'{key}_mrnas'], counts[f'{key}_loci'] = int(m.group(2)), int(m.group(3))
                continue
            m = GFFCOMPARE_MATCH_RE.match(line)
            if m:
                counts['matching_' + m.group(1).replace(' ', '_')] = int(m.group(2))
                continue
            m = GFFCOMPARE_MISSED_RE.match(line)
            if m:
                key = f'{m.group(1).lower()}_{m.group(2)}'
                counts[key], counts[f'{key}_total'] = int(m.group(3)), int(m.group(4))
    levels = {level: {'sn': sn, 'pr': pr} for level, (sn, pr) in parse_accuracy(stats_path).items()}
    return {'levels': levels, 'counts': counts}


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
