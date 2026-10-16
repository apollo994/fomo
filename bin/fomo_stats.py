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
import re
import sys
from pathlib import Path

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

# Pseudo-source ids of aggregate models (every source pooled), which per-source tables skip.
# `from_<id>` filename convention, so per-source tables must skip them.
AGGREGATE_IDS = ("allModels_raw",)


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


def f1(sn: float, pr: float) -> float:
    return 0.0 if (sn + pr) == 0 else 2.0 * sn * pr / (sn + pr)


# Curation report (curate_models.py), its data record and the tables rendered from it.
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


# ── gff-feature-stats rows ────────────────────────────────────────────────────────
# The table rows of one gff-feature-stats JSON, formatted for a TSV. Blocks the tool omits are NA,
# not 0: intron_stats is absent when every transcript is single-exon, cds_stats for non-coding.
NA = "NA"

# Column order for the per-transcript-type table.
GFFSTATS_TRANSCRIPT_COLUMNS: list[str] = [
    "role",
    "kind",
    "transcript_type",
    "n_genes",
    "n_transcripts",
    "n_exons",
    "mean_exon_length",
    "max_exon_length",
    "n_introns",
    "mean_intron_length",
    "min_intron_length",
    "max_intron_length",
    "mean_transcript_length",
    "min_transcript_length",
    "max_transcript_length",
    "mean_spliced_length",
    "min_spliced_length",
    "max_spliced_length",
    "n_cds",
    "mean_cds_length",
    "top_biotype",
]

# Column order for the per-gene-category table.
GFFSTATS_GENE_COLUMNS: list[str] = [
    "role",
    "kind",
    "gene_category",
    "n_genes",
    "mean_gene_length",
    "min_gene_length",
    "max_gene_length",
    "n_transcript_types",
    "top_biotype",
]


def fmt(v) -> str:
    """Render a metric for the TSV. None → NA; integral floats lose the .0."""
    if v is None:
        return NA
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        if v != v:  # NaN
            return NA
        return str(int(v)) if v.is_integer() else f"{v:.2f}"
    return str(v)


def length_field(block: dict | None, key: str) -> object:
    """Pull min/max/mean out of a `{min, max, mean}` sub-object, or None."""
    if not isinstance(block, dict):
        return None
    value = block.get(key)
    return value if isinstance(value, int | float) else None


def top_key(counts) -> str | None:
    """Highest-count key of a `{name: count}` map, ties broken by name."""
    if not isinstance(counts, dict) or not counts:
        return None
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def gffstats_transcript_rows(stats: dict, role: str, kind: str) -> list[dict]:
    """Build one row per transcript type, in the JSON's own (count-desc) order."""
    rows: list[dict] = []
    for ttype, ts in (stats.get("transcript_type_stats") or {}).items():
        length_stats = ts.get("length_stats")
        exon = ts.get("exon_stats") or {}
        exon_len = exon.get("length")
        spliced = exon.get("concatenated_length")
        intron = ts.get("intron_stats") or {}
        intron_len = intron.get("length")
        cds = ts.get("cds_stats") or {}
        cds_len = cds.get("length")
        assoc = ts.get("associated_genes") or {}

        rows.append(
            {
                "role": role,
                "kind": kind,
                "transcript_type": ttype,
                "n_genes": fmt(assoc.get("total_count")),
                "n_transcripts": fmt(ts.get("total_count")),
                "n_exons": fmt(exon.get("total_count")),
                "mean_exon_length": fmt(length_field(exon_len, "mean")),
                "max_exon_length": fmt(length_field(exon_len, "max")),
                "n_introns": fmt(intron.get("total_count")),
                "mean_intron_length": fmt(length_field(intron_len, "mean")),
                "min_intron_length": fmt(length_field(intron_len, "min")),
                "max_intron_length": fmt(length_field(intron_len, "max")),
                "mean_transcript_length": fmt(length_field(length_stats, "mean")),
                "min_transcript_length": fmt(length_field(length_stats, "min")),
                "max_transcript_length": fmt(length_field(length_stats, "max")),
                "mean_spliced_length": fmt(length_field(spliced, "mean")),
                "min_spliced_length": fmt(length_field(spliced, "min")),
                "max_spliced_length": fmt(length_field(spliced, "max")),
                "n_cds": fmt(cds.get("total_count")),
                "mean_cds_length": fmt(length_field(cds_len, "mean")),
                "top_biotype": fmt(top_key(ts.get("biotype_counts"))),
            }
        )

    if not rows:
        placeholder = {c: NA for c in GFFSTATS_TRANSCRIPT_COLUMNS}
        placeholder.update({"role": role, "kind": kind, "transcript_type": "none"})
        rows.append(placeholder)
    return rows


def gffstats_gene_rows(stats: dict, role: str, kind: str) -> list[dict]:
    """Build one row per gene category (coding / pseudogene / non_coding)."""
    rows: list[dict] = []
    for category, gc in (stats.get("gene_category_stats") or {}).items():
        length_stats = gc.get("length_stats")
        ttype_counts = gc.get("transcript_type_counts")
        rows.append(
            {
                "role": role,
                "kind": kind,
                "gene_category": category,
                "n_genes": fmt(gc.get("total_count")),
                "mean_gene_length": fmt(length_field(length_stats, "mean")),
                "min_gene_length": fmt(length_field(length_stats, "min")),
                "max_gene_length": fmt(length_field(length_stats, "max")),
                "n_transcript_types": fmt(len(ttype_counts) if isinstance(ttype_counts, dict) else None),
                "top_biotype": fmt(top_key(gc.get("biotype_counts"))),
            }
        )

    if not rows:
        placeholder = {c: NA for c in GFFSTATS_GENE_COLUMNS}
        placeholder.update({"role": role, "kind": kind, "gene_category": "none"})
        rows.append(placeholder)
    return rows
