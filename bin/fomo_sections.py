#!/usr/bin/env python3
"""
The MultiQC sections of the per-target report: for each, the records that feed it, its rows
and their labels, and the custom-content file name its `sp` pattern matches
(assets/multiqc/sections.yml; tests/test_multiqc_registry.py keeps the two consistent).

Row labels turn a record's identity into a display name. The order tokens in them
(`1_raw`, `2_lncRNA`, ...) sort the rows of a table and exist only here.

A module, not a script; render_mqc.py uses it.
"""

import colorsys
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fomo_stats  # noqa: E402


def cell(value) -> str:
    """A table cell: None → NA, booleans yes/no, numbers as written."""
    if value is None:
        return fomo_stats.NA
    if isinstance(value, bool):
        return 'yes' if value else 'no'
    return value if isinstance(value, str) else repr(value) if isinstance(value, float) else str(value)


def annotation_label(ident: dict) -> str:
    """<species>.<order token>, with `_target` for a target's own annotation."""
    token = '1_raw' if ident['kind'] == 'raw' else f'{ident["order"] + 1}_{ident["label"]}'
    return f'{ident["species"]}.{token}' + ('_target' if ident['side'] == 'target' else '')


def track_label(ident: dict) -> str:
    return f'{ident["order"]}_{ident["label"]}'


# ── row builders: record → (file name stem, [row]) with Sample first ────────────
def gffstats_input(rec: dict):
    ident, base = rec['identity'], annotation_label(rec['identity'])
    rows = fomo_stats.gffstats_transcript_rows(rec['data'], ident['side'], ident['kind'])
    return base, [{'Sample': f'{base}.{r["transcript_type"]}', **r} for r in rows]


def gffstats_genes(rec: dict):
    ident, base = rec['identity'], annotation_label(rec['identity'])
    rows = fomo_stats.gffstats_gene_rows(rec['data'], ident['side'], ident['kind'])
    return base, [{'Sample': f'{base}.{r["gene_category"]}', **r} for r in rows]


def gffstats_projection(rec: dict):
    ident = rec['identity']
    base = f'{ident["target"]}.from_{ident["source"]}.{track_label(ident)}'
    rows = fomo_stats.gffstats_transcript_rows(rec['data'], 'projection', 'projected')
    return base, [{'Sample': f'{base}.{r["transcript_type"]}', **r} for r in rows]


def fasta_stats(rec: dict):
    ident = rec['identity']
    base = f'{ident["species"]}.{ident["track"]}'
    return base, [{'Sample': base} | {k: cell(v) for k, v in rec['data'].items()}]


ALIGNMENT_COLUMNS = (
    'raw_total_sequences',
    'reads_mapped',
    'reads_mapped_percent',
    'non_primary_alignments',
    'error_rate',
    'average_quality',
    'average_length',
)


def samtools_align(rec: dict):
    ident = rec['identity']
    base = f'{ident["target"]}.allModels.{track_label(ident)}'
    return base, [{'Sample': base} | {k: cell(rec['data'].get(k)) for k in ALIGNMENT_COLUMNS}]


def transcript_filter(rec: dict):
    ident, data = rec['identity'], rec['data']
    base = f'{ident["species"]}.{ident["track"]}'
    counts = {k: cell(data[k]) for k in ('n_transcripts', 'n_stranded', 'n_spliced', 'n_selected')}
    return base, [{'Sample': base, 'feature': ident['track']} | counts]


def td2_coding(rec: dict):
    ident, data = rec['identity'], rec['data']
    base = ident['species'] if ident['stage'] == 'input' else f'{ident["target"]}.allModels'
    pct = 100.0 * data['n_coding'] / data['n_in'] if data['n_in'] else 0.0
    row = {'Sample': base, 'stage': ident['stage'], 'feature': ident['label']}
    row |= {k: cell(data[k]) for k in ('n_in', 'n_coding', 'n_kept')} | {'pct_coding': f'{pct:.2f}'}
    return f'{base}.{ident["stage"]}', [row]


def curation_base(rec: dict) -> str:
    return f'{rec["identity"]["target"]}.{rec["identity"]["track"]}'


def curation_observed(rec: dict):
    observed = rec['data']['observed']
    base = curation_base(rec)
    return base, [{'Sample': base} | {k: cell(observed.get(k, 0)) for k in fomo_stats.CURATION_OBSERVED}]


def curation_discarded(rec: dict):
    data, base = rec['data'], curation_base(rec)
    has_ref, discarded = data['has_reference'], data['discarded_by_reference']
    row = {'Sample': base, 'reference': 'yes' if has_ref else 'no'}
    row |= {k: cell(discarded.get(k, 0) if has_ref else None) for k in fomo_stats.CURATION_DISCARDED}
    prefix = 'chains_lost_by_primary_type:'
    row |= {f'chains_lost:{k[len(prefix) :]}': cell(v) for k, v in discarded.items() if k.startswith(prefix)}
    return base, [row]


def curation_novel(rec: dict):
    data, base = rec['data'], curation_base(rec)
    novel = data['novel']
    row = {'Sample': base} | {k: cell(novel.get(k, 0)) for k in fomo_stats.CURATION_NOVEL}
    row |= {f'location:{c}': cell(novel.get(f'genes_location:{c}', 0)) for c in fomo_stats.CURATION_LOCATION_COLUMNS}
    row |= {f'median_{m}': cell(novel.get(f'{m}_median')) for m in fomo_stats.CURATION_MEDIANS}
    distance = novel.get('intergenic_distance_median') if data['has_reference'] else None
    return base, [row | {'median_intergenic_distance': cell(distance)}]


@dataclass(frozen=True)
class Section:
    id: str  # MultiQC section id: the custom_data and sp key
    record: str  # record section that feeds it
    select: Callable[[dict], bool]  # which records, by identity
    rows: Callable[[dict], tuple]  # record → (file stem, rows)

    @property
    def suffix(self) -> str:
        return f'_{self.id}_mqc.tsv'


def always(_ident: dict) -> bool:
    return True


SECTIONS = (
    Section('gffstats_input', 'gffstats', lambda i: i['kind'] != 'projected', gffstats_input),
    Section('gffstats_genes', 'gffstats', lambda i: i['kind'] != 'projected', gffstats_genes),
    Section('gffstats_projection', 'gffstats', lambda i: i['kind'] == 'projected', gffstats_projection),
    Section('seqkit_stats', 'fasta_stats', always, fasta_stats),
    Section('samtools_align', 'alignment', always, samtools_align),
    Section('transcript_filter', 'transcript_filter', always, transcript_filter),
    Section('td2_coding', 'td2_coding', always, td2_coding),
    Section('curation_observed', 'curation', always, curation_observed),
    Section('curation_discarded', 'curation', always, curation_discarded),
    Section('curation_novel', 'curation', always, curation_novel),
)

# ── accuracy scatter: Sn vs Pr per level, from the accuracy records ─────────────
SCATTER_ID = 'gffcompare_accuracy_custom'
SCATTER_SUFFIX = '.gffcompare_accuracy_mqc.json'
SYMBOL_BY_CLASS = {'lncRNA': 'circle', 'mRNA': 'square', 'decoy': 'triangle-up'}
AGGREGATE_COLOUR = '#222222'


def source_colours(sources) -> dict:
    """The aggregates in one fixed dark colour; every species its own hue, evenly spaced, so the
    colours never repeat however many sources a run has."""
    species = sorted(s for s in sources if s not in fomo_stats.AGGREGATE_IDS)
    colours = {s: AGGREGATE_COLOUR for s in sources if s in fomo_stats.AGGREGATE_IDS}
    for i, source in enumerate(species):
        r, g, b = colorsys.hls_to_rgb(i / max(len(species), 1), 0.45, 0.65)
        colours[source] = f'#{round(r * 255):02x}{round(g * 255):02x}{round(b * 255):02x}'
    return colours


def scatter(records: list[dict]) -> dict:
    colours = source_colours({rec['identity']['source'] for rec in records})
    data = []
    for level in fomo_stats.LEVELS:
        points = {}
        for rec in records:
            ident, cell_ = rec['identity'], rec['data']['levels'].get(level)
            if cell_ is None:
                continue
            name = f'{ident["source"]}.{ident["track"]}'
            klass = 'decoy' if ident['track'].endswith('.decoy') else ident['label']
            points[name] = {
                'x': round(cell_['sn'] / 100.0, 4),
                'y': round(cell_['pr'] / 100.0, 4),
                'color': colours[ident['source']],
                'marker_symbol': SYMBOL_BY_CLASS[klass],
                'name': name,
            }
        data.append(points)
    return {
        'id': SCATTER_ID,
        'section_name': 'Accuracy values',
        'description': (
            'Sensitivity vs Precision per feature level (use the dropdown to switch level). Marker symbol = '
            'feature class (lncRNA ●, mRNA ■, decoy ▲); colour = source species, or dark grey for allModels '
            '(every source pooled, raw). Hover a point for its name.'
        ),
        'plot_type': 'scatter',
        'pconfig': {
            'id': f'{SCATTER_ID}_plot',
            'title': 'GffCompare: Accuracy values',
            'xlab': 'Sensitivity',
            'ylab': 'Precision',
            'xmin': 0,
            'xmax': 1,
            'ymin': 0,
            'ymax': 1,
            'data_labels': [
                {'name': level.replace('_', ' '), 'xlab': 'Sensitivity', 'ylab': 'Precision'}
                for level in fomo_stats.LEVELS
            ],
        },
        'data': data,
    }
