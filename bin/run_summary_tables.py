#!/usr/bin/env python3
"""
Digest a whole FOMO run into run-level MultiQC sections.

The per-target reports answer "how did this target do". This answers "what did this run do":
species per role, how much annotation went in and how much survived filtering, how much each
target gained, and which donors are worth using.

Everything comes from the run's data records (fomo-data/1, one directory) and the roles CSV, so
this script adds no measurement, only cross-target arithmetic; identity comes from each record's
identity, never from a file name or a row label. Nothing here parses a GFF or a BAM.

Outputs, in report order (see assets/multiqc/run_summary_main.yml):
  run_overview_mqc.tsv            run overview: roles, pair counts, params
  run_species_mqc.tsv             per species: role + its own annotation input
  run_source_funnel_mqc.{tsv,json} per source: raw → stranded → spliced → longest → non-coding funnel
  run_target_types_mqc.tsv        per (target × transcript type): genes before / curated / after
  run_target_before_after_mqc.json  the lncRNA genes, as a stacked bargraph
  run_accuracy_heatmap_{f1,sn,pr}_mqc.json  source × target accuracy matrices
  run_aggregate_accuracy_mqc.{tsv,json}  allModels (every source, raw) vs best single donor
  run_donor_ranking_mqc.tsv       per source: F1 across targets
  run_summary.json                every number above, machine-readable
"""

import argparse
import csv
import json
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fomo_io  # noqa: E402
import fomo_stats  # noqa: E402
from fomo_stats import AGGREGATE_IDS, f1, fmt  # noqa: E402

log = fomo_io.get_logger('run_summary')
GENE_CATEGORIES = ('coding', 'pseudogene', 'non_coding')


# ── small helpers ───────────────────────────────────────────────────────────────
def num(value) -> float | None:
    """A table cell or record value as a float, or None for NA/missing/unparseable."""
    if value is None or value in (fomo_stats.NA, '', 'None'):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def pct(part: float | None, whole: float | None) -> float | None:
    if part is None or whole in (None, 0):
        return None
    return 100.0 * part / whole


def add(*values: float | None) -> float | None:
    """Sum that propagates None: 'unknown + 5' is unknown, not 5."""
    if any(v is None for v in values):
        return None
    return sum(values)


def track_of(feature_type: str, decoy: bool) -> str:
    return feature_type + ('.decoy' if decoy else '')


def write_tsv(name: str, header, rows) -> None:
    with open(name, 'w', encoding='utf-8') as fh:
        fh.write('\t'.join(header) + '\n')
        for row in rows:
            fh.write('\t'.join(fmt(cell) for cell in row) + '\n')
    log.info(f'wrote {name} ({len(rows)} rows)')


def write_json(name: str, doc: dict) -> None:
    with open(name, 'w', encoding='utf-8') as fh:
        json.dump(doc, fh, indent=2)
    log.info(f'wrote {name}')


# ── inputs ──────────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class RunContext:
    """The run's species roles and switches, derived once."""

    roles: dict
    feature_types: tuple
    include_decoy: bool
    psauron_min: float

    @property
    def sources(self) -> list[str]:
        return sorted(s for s, r in self.roles.items() if r['is_source'])

    @property
    def targets(self) -> list[str]:
        return sorted(s for s, r in self.roles.items() if r['is_target'])


def read_roles(path: Path) -> dict:
    roles = {}
    with path.open(encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            role = row['role']
            roles[row['species']] = {
                'role': role,
                'has_gff3': row['has_gff3'] == 'yes',
                'is_source': role in ('source', 'both'),
                'is_target': role in ('target', 'both'),
            }
    if not roles:
        fomo_io.die(f'no species rows in {path}')
    return roles


class RecordIndex:
    """The run's records, indexed by identity. Every index is built once, on first use."""

    def __init__(self, directory: Path):
        self.by_section = defaultdict(list)
        for record in fomo_stats.read_records(directory):
            self.by_section[record['section']].append(record)
        log.info('staged ' + ', '.join(f'{len(v)} {k}' for k, v in sorted(self.by_section.items())) + ' records')

    @cached_property
    def _gffstats(self) -> dict:
        """(species, side, kind, track) -> record, for the annotation (not projected) stats."""
        index = {}
        for rec in self.by_section['gffstats']:
            ident = rec['identity']
            if ident['kind'] != 'projected':
                index[(ident['species'], ident['side'], ident['kind'], ident.get('track'))] = rec
        return index

    def _annotation_record(self, species: str, side: str, kind: str, track: str | None):
        """A target's raw annotation may only have been measured on the source side."""
        rec = self._gffstats.get((species, side, kind, track))
        if rec is None and side == 'target' and kind == 'raw':
            rec = self._gffstats.get((species, 'source', kind, track))
        return rec

    def transcript_rows(self, species: str, side: str, kind: str, track: str | None = None) -> dict:
        """{transcript type: table row} of an annotation's gff-feature-stats."""
        rec = self._annotation_record(species, side, kind, track)
        if rec is None:
            return {}
        rows = fomo_stats.gffstats_transcript_rows(rec['data'], side, kind)
        return {row['transcript_type']: row for row in rows}

    def gene_rows(self, species: str, side: str, kind: str) -> dict:
        """{gene category: table row} of an annotation's gff-feature-stats."""
        rec = self._annotation_record(species, side, kind, None)
        if rec is None:
            return {}
        return {row['gene_category']: row for row in fomo_stats.gffstats_gene_rows(rec['data'], side, kind)}

    @cached_property
    def projected_models(self) -> dict:
        """(target, model, track) -> projected transcripts (every transcript type)."""
        index = {}
        for rec in self.by_section['gffstats']:
            ident = rec['identity']
            if ident['kind'] == 'projected':
                rows = fomo_stats.gffstats_transcript_rows(rec['data'], 'projection', 'projected')
                index[(ident['target'], ident['source'], ident['track'])] = add(
                    *[num(r.get('n_transcripts')) or 0.0 for r in rows]
                )
        return index

    @cached_property
    def td2(self) -> dict:
        """(stage, species) -> counts; the species is the target for the projected stage."""
        return {(r['identity']['stage'], r['identity']['species']): r['data'] for r in self.by_section['td2_coding']}

    @cached_property
    def transcript_filter(self) -> dict:
        """(species, track) -> FILTER_TRANSCRIPT funnel."""
        return {
            (r['identity']['species'], r['identity']['track']): r['data'] for r in self.by_section['transcript_filter']
        }

    @cached_property
    def alignment(self) -> dict:
        """(target, track) -> samtools stats."""
        return {(r['identity']['target'], r['identity']['track']): r['data'] for r in self.by_section['alignment']}

    @cached_property
    def accuracy(self) -> dict:
        """(target, source, track) -> {sn, pr, f1} at transcript level."""
        acc = {}
        for rec in self.by_section['accuracy']:
            ident, cell = rec['identity'], rec['data']['levels'].get('Transcript')
            if cell is None:
                log.warning(f'no Transcript level in the accuracy record of {ident}, skipping')
                continue
            acc[(ident['target'], ident['source'], ident['track'])] = {
                'sn': cell['sn'],
                'pr': cell['pr'],
                'f1': f1(cell['sn'], cell['pr']),
            }
        return acc

    @cached_property
    def curation(self) -> dict:
        """(target, track) -> curation report (observed / discarded / novel)."""
        return {(r['identity']['target'], r['identity']['track']): r['data'] for r in self.by_section['curation']}


# ── sections: compute (pure) and emit ──────────────────────────────────────────
def compute_overview(ctx: RunContext) -> tuple[list, dict]:
    roles = ctx.roles
    sources, targets = ctx.sources, ctx.targets
    # Benchmarked = donor AND target with a gff3 (role both): the gffcompare reference is the
    # species' filtered annotation, which only donors have. A pure target may carry a gff3 as its
    # curation reference only.
    annotated = [t for t in targets if roles[t]['has_gff3'] and roles[t]['is_source']]
    with_reference = [t for t in targets if roles[t]['has_gff3']]
    # A species is never projected onto itself.
    n_pairs = len(sources) * len(targets) - sum(1 for s in sources if s in targets)
    n_role = {role: sum(1 for r in roles.values() if r['role'] == role) for role in ('source', 'target', 'both')}
    feature_types = list(ctx.feature_types)
    metrics = [
        ('Species in samplesheet', len(roles)),
        ('… role source (donor only)', n_role['source']),
        ('… role target (projected onto, not a donor)', n_role['target']),
        ('… role both (donor and target)', n_role['both']),
        ('Sources (S)', len(sources)),
        ('Targets (T)', len(targets)),
        ('Targets with a reference annotation (curation)', len(with_reference)),
        ('Targets benchmarked (role both)', len(annotated)),
        ('Source × target pairs (no self-projection)', n_pairs),
        ('Feature types transferred', ', '.join(feature_types)),
        ('Decoy track (--include_decoy)', ctx.include_decoy),
        ('PSAURON cutoff (--td2_psauron_min)', ctx.psauron_min),
    ]
    overview = {
        'n_species': len(roles),
        'n_source_only': n_role['source'],
        'n_target_only': n_role['target'],
        'n_both': n_role['both'],
        'n_sources': len(sources),
        'n_targets': len(targets),
        'n_targets_annotated': len(annotated),
        'n_targets_with_reference': len(with_reference),
        'n_pairs': n_pairs,
        'feature_types': feature_types,
        'include_mrna': 'mRNA' in feature_types,
        'include_decoy': ctx.include_decoy,
        'td2_psauron_min': ctx.psauron_min,
    }
    return metrics, overview


SPECIES_HEADER = [
    'Sample',
    'role',
    'has_gff3',
    'is_source',
    'is_target',
    'benchmarked',
    'input_genes_total',
    'input_lncRNA_transcripts_raw',
    'spliced_kept',
    'td2_coding_dropped',
    'lncRNA_donated',
]


def own_annotation(idx: RecordIndex, species: str, info: dict) -> tuple[dict, float | None]:
    """A species' own raw lncRNA row and gene total, from the source side if it donates, else
    the target side (a both species has both, and they are the same annotation)."""
    side = 'source' if info['is_source'] else 'target'
    genes = idx.gene_rows(species, side, 'raw')
    total = None
    for category in GENE_CATEGORIES:
        if category in genes:
            total = add(total or 0.0, num(genes[category].get('n_genes')))
    return idx.transcript_rows(species, side, 'raw').get('lnc_RNA'), total


def compute_species(ctx: RunContext, idx: RecordIndex) -> dict:
    records = {}
    for species in sorted(ctx.roles):
        info = ctx.roles[species]
        raw_row, genes_total = own_annotation(idx, species, info)
        td2 = idx.td2.get(('input', species)) or {}
        records[species] = {
            'role': info['role'],
            'has_gff3': info['has_gff3'],
            'is_source': info['is_source'],
            'is_target': info['is_target'],
            'benchmarked': info['is_target'] and info['is_source'] and info['has_gff3'],
            'input_genes_total': genes_total,
            'input_lncRNA_transcripts_raw': num((raw_row or {}).get('n_transcripts')),
            'spliced_kept': num(td2.get('n_in')),
            'td2_coding_dropped': num(td2.get('n_coding')),
            'lncRNA_donated': num(td2.get('n_kept')),
        }
    return records


FUNNEL_HEADER = [
    'Sample',
    'raw_lncRNA_transcripts',
    'dropped_unstranded',
    'dropped_single_exon',
    'dropped_non_longest',
    'dropped_coding',
    'kept',
    'pct_kept',
    'n_genes_raw',
    'n_genes_kept',
    'mean_transcript_length_kept',
    'mean_spliced_length_kept',
    'mean_exon_length_kept',
    'exons_per_transcript_kept',
]


def difference(a: float | None, b: float | None) -> float | None:
    return None if a is None or b is None else a - b


def check_funnel_consistency(species: str, raw, n_tx, selected, td2_in, kept, gff_kept) -> None:
    """Each step must hand the next exactly what it says it kept; a mismatch means two artefacts
    stopped describing the same set."""
    checks = (
        (selected, td2_in, f'FILTER_TRANSCRIPT selected {selected} lncRNA but TD2 received {td2_in}'),
        (n_tx, raw, f'FILTER_TRANSCRIPT saw {n_tx} lncRNA but the raw GFF stats count {raw}'),
        (kept, gff_kept, f'TD2 kept {kept} lncRNA but the filtered GFF3 has {gff_kept}'),
    )
    for a, b, message in checks:
        if a is not None and b is not None and abs(a - b) > 0.5:
            log.warning(f'{species}: {message}')


def compute_source_funnel(ctx: RunContext, idx: RecordIndex) -> tuple[dict, dict]:
    records, plot = {}, {}
    for species in ctx.sources:
        raw_row = idx.transcript_rows(species, 'source', 'raw').get('lnc_RNA') or {}
        kept_row = idx.transcript_rows(species, 'source', 'filtered', 'lnc_RNA').get('lnc_RNA') or {}
        td2 = idx.td2.get(('input', species)) or {}
        funnel = idx.transcript_filter.get((species, 'lnc_RNA')) or {}
        raw, coding, kept = num(raw_row.get('n_transcripts')), num(td2.get('n_coding')), num(td2.get('n_kept'))
        n_tx, stranded = num(funnel.get('n_transcripts')), num(funnel.get('n_stranded'))
        spliced, selected = num(funnel.get('n_spliced')), num(funnel.get('n_selected'))
        unstranded, single_exon = difference(n_tx, stranded), difference(stranded, spliced)
        non_longest = difference(spliced, selected)
        gff_kept = num(kept_row.get('n_transcripts'))
        check_funnel_consistency(species, raw, n_tx, selected, num(td2.get('n_in')), kept, gff_kept)
        n_exons = num(kept_row.get('n_exons'))
        records[species] = {
            'raw_lncRNA_transcripts': raw,
            'dropped_unstranded': unstranded,
            'dropped_single_exon': single_exon,
            'dropped_non_longest': non_longest,
            'dropped_coding': coding,
            'kept': kept,
            'pct_kept': pct(kept, raw),
            'n_genes_raw': num(raw_row.get('n_genes')),
            'n_genes_kept': num(kept_row.get('n_genes')),
            'mean_transcript_length_kept': num(kept_row.get('mean_transcript_length')),
            'mean_spliced_length_kept': num(kept_row.get('mean_spliced_length')),
            'mean_exon_length_kept': num(kept_row.get('mean_exon_length')),
            'exons_per_transcript_kept': None if (n_exons is None or not gff_kept) else n_exons / gff_kept,
        }
        # Stacked bars only make sense when the whole funnel resolved.
        if None not in (kept, unstranded, single_exon, non_longest, coding):
            plot[species] = {
                'kept': kept,
                'dropped_unstranded': unstranded,
                'dropped_single_exon': single_exon,
                'dropped_non_longest': non_longest,
                'dropped_coding': coding,
            }
    return records, plot


def emit_source_funnel(records: dict, plot: dict) -> None:
    write_tsv(
        'run_source_funnel_mqc.tsv',
        FUNNEL_HEADER,
        [[s] + [r[k] for k in FUNNEL_HEADER[1:]] for s, r in records.items()],
    )
    write_json(
        'run_source_funnel_mqc.json',
        {
            'id': 'run_source_funnel_plot',
            'section_name': 'Source lncRNA — what survived filtering',
            'description': (
                'Per source species: lncRNA transcripts donated (kept) and the four '
                'reasons the rest were dropped — no consistent strand ("?" or '
                'trans-spliced across both strands; almost always zero), single-exon (FOMO transfers spliced '
                'models only), not the longest spliced isoform of its gene (one '
                'transcript per gene; zero with --longest_isoform false) and coding '
                'potential (a complete ORF with PSAURON score at or above the cutoff). '
                'Bars sum to the raw lncRNA count.'
            ),
            'plot_type': 'bargraph',
            'pconfig': {
                'id': 'run_source_funnel_bargraph',
                'title': 'FOMO: source lncRNA filtering funnel',
                'ylab': 'transcripts',
                'cpswitch_counts_label': 'Transcripts',
            },
            'data': plot,
        },
    )


def target_header(ctx: RunContext) -> list[str]:
    header = ['Sample', 'target', 'transcript_type', 'genes_before', 'curated_genes', 'genes_after', 'pct_increase']
    header.append('reads_mapped_percent')
    if ctx.include_decoy:
        # A decoy curated gene is a chance cross-species agreement; normalised per multi-exon
        # input model, because decoy_cap makes the decoy and real inputs different sizes.
        header += ['curated_decoy_genes', 'est_curated_fdr_pct']
    return header


def curated(idx: RecordIndex, target: str, track: str, table: str = 'novel', key: str = 'genes') -> float | None:
    report = idx.curation.get((target, track))
    return None if report is None else num(report[table].get(key, 0))


def decoy_columns(idx: RecordIndex, target: str, ttype: str, added: float | None, enabled: bool) -> dict:
    decoy = curated(idx, target, track_of(ttype, True)) if enabled else None
    real_in = curated(idx, target, ttype, 'observed', 'multi_exon_models')
    decoy_in = curated(idx, target, track_of(ttype, True), 'observed', 'multi_exon_models')
    rate_real = None if added is None or not real_in else added / real_in
    rate_decoy = None if decoy is None or not decoy_in else decoy / decoy_in
    fdr = None if rate_real in (None, 0) or rate_decoy is None else 100.0 * rate_decoy / rate_real
    return {'curated_decoy_genes': decoy, 'est_curated_fdr_pct': fdr}


def compute_targets(ctx: RunContext, idx: RecordIndex) -> dict:
    """Per (target × transcript type): genes already annotated, curated genes FOMO adds and the
    sum, counted in genes on both sides. Only the enabled tracks are curated: a type FOMO does not
    transfer gains a known zero, an enabled type without a curation report is unknown (NA)."""
    records = {}
    for target in ctx.targets:
        records[target] = {}  # every target appears, even one without any transcript type
        has_ref = ctx.roles[target]['has_gff3']
        raw = idx.transcript_rows(target, 'target', 'raw')
        types = set(raw) | {track for (t, track) in idx.curation if t == target and not track.endswith('.decoy')}
        for ttype in sorted(types):
            enabled = ttype in ctx.feature_types
            before = num((raw.get(ttype) or {}).get('n_genes')) if has_ref else None
            added = curated(idx, target, ttype) if enabled else 0.0
            mapped = num((idx.alignment.get((target, ttype)) or {}).get('reads_mapped_percent')) if enabled else None
            record = {
                'target': target,
                'transcript_type': ttype,
                'genes_before': before,
                'curated_genes': added,
                'genes_after': add(before, added) if has_ref else added,
                'pct_increase': pct(added, before),
                'reads_mapped_percent': mapped,
            }
            if ctx.include_decoy:
                record |= decoy_columns(idx, target, ttype, added, enabled)
            records[target][ttype] = record
    return records


def emit_targets(ctx: RunContext, records: dict) -> None:
    header = target_header(ctx)
    rows = [[f'{t}.{tt}'] + [r[k] for k in header[1:]] for t, by_type in records.items() for tt, r in by_type.items()]
    write_tsv('run_target_types_mqc.tsv', header, rows)
    bar = {}
    for target, by_type in records.items():
        lnc = by_type.get('lnc_RNA', {})
        if lnc.get('curated_genes') is not None:
            bar[target] = {
                'existing_lncRNA_genes': lnc.get('genes_before') or 0.0,
                'curated_lncRNA_genes': lnc['curated_genes'],
            }
    write_json(
        'run_target_before_after_mqc.json',
        {
            'id': 'run_target_before_after',
            'section_name': 'Targets — lncRNA genes before and after FOMO',
            'description': (
                'Per target: lncRNA genes already in its reference annotation versus the curated '
                'lncRNA genes FOMO adds — intron chains shared exactly by at least '
                '--curate_min_species source species, none overlapping a reference exon. A '
                'target without a reference annotation starts at zero.'
            ),
            'plot_type': 'bargraph',
            'pconfig': {
                'id': 'run_target_before_after_bargraph',
                'title': 'FOMO: target lncRNA genes before and after curation',
                'ylab': 'genes',
            },
            'data': bar,
        },
    )


# One section per metric, not one plot with a dataset switcher: MultiQC 1.35 accepts
# multi-dataset bargraphs but not heatmaps.
HEATMAP_METRICS = [
    ('f1', 'F1', 'F1 (the harmonic mean of the two below)'),
    ('sn', 'Sensitivity', "Sensitivity — how much of the target's real annotation the projection recovered"),
    ('pr', 'Precision', 'Precision — how much of the projection is real annotation rather than noise'),
]
HEATMAP_NOTE = (
    'lncRNA transcript-level gffcompare accuracy for each source → target '
    'projection. A blank cell means that pair was not benchmarked — the target '
    'is not benchmarked (only role `both` targets are). The diagonal is always '
    'blank: a species is never projected onto itself (plans/23).'
)


def compute_heatmap(ctx: RunContext, idx: RecordIndex) -> dict:
    """metric -> source × target matrix (None where not benchmarked)."""
    matrices = {}
    for metric, _, _ in HEATMAP_METRICS:
        matrices[metric] = [
            [
                None if cell is None else round(cell[metric], 2)
                for cell in (idx.accuracy.get((target, source, 'lnc_RNA')) for target in ctx.targets)
            ]
            for source in ctx.sources
        ]
    return matrices


def emit_heatmap(ctx: RunContext, matrices: dict) -> None:
    for metric, label, blurb in HEATMAP_METRICS:
        write_json(
            f'run_accuracy_heatmap_{metric}_mqc.json',
            {
                'id': f'run_accuracy_heatmap_{metric}',
                'section_name': f'Accuracy {label} — every source against every target',
                'description': f'{blurb}. {HEATMAP_NOTE}',
                'plot_type': 'heatmap',
                'xcats': ctx.targets,
                'ycats': ctx.sources,
                'pconfig': {
                    'id': f'run_accuracy_heatmap_{metric}_plot',
                    'title': f'FOMO: source × target lncRNA {label}',
                    'xlab': 'target',
                    'ylab': 'source',
                    'min': 0,
                    'max': 100,
                },
                'data': matrices[metric],
            },
        )


def aggregate_header(ctx: RunContext) -> list[str]:
    header = ['Sample', 'target', 'feature_type', 'model', 'best_source', 'sn', 'pr', 'f1', 'n_models']
    header.append('delta_f1_vs_best_single')
    return header + (['decoy_models', 'decoy_f1', 'est_fdr_pct'] if ctx.include_decoy else [])


def best_single_source(idx: RecordIndex, target: str, ft: str):
    candidates = {
        source: cell['f1']
        for (tgt, source, track), cell in idx.accuracy.items()
        if tgt == target and track == ft and source not in AGGREGATE_IDS and source != target
    }
    best = max(candidates, key=lambda s: (candidates[s], s), default=None)
    return best, candidates.get(best) if best else None


def aggregate_record(
    ctx: RunContext, idx: RecordIndex, target: str, ft: str, label: str, model: str, best
) -> dict | None:
    best_source, best_f1 = best
    cell = idx.accuracy.get((target, model, ft))
    if cell is None:
        return None
    real_models = idx.projected_models.get((target, model, ft))
    record = {
        'target': target,
        'feature_type': ft,
        'model': label,
        'best_source': best_source,
        'sn': cell['sn'],
        'pr': cell['pr'],
        'f1': cell['f1'],
        'n_models': real_models,
        'delta_f1_vs_best_single': None if best_f1 is None else cell['f1'] - best_f1,
    }
    if ctx.include_decoy:
        decoy_models = idx.projected_models.get((target, model, track_of(ft, True)))
        decoy_cell = idx.accuracy.get((target, model, track_of(ft, True)))
        record['decoy_models'] = decoy_models
        record['decoy_f1'] = None if decoy_cell is None else decoy_cell['f1']
        record['est_fdr_pct'] = pct(decoy_models, real_models)
    return record


def compute_aggregates(ctx: RunContext, idx: RecordIndex) -> dict:
    """Per (target, feature type): the best real donor and each pooled model, side by side. The
    row label is a category name in the bargraph, so it does not carry the species."""
    records = defaultdict(dict)
    for target in ctx.targets:
        for ft in ctx.feature_types:
            best = best_single_source(idx, target, ft)
            entries = ([('best_single_source', best[0])] if best[0] else []) + [(m, m) for m in AGGREGATE_IDS]
            for label, model in entries:
                record = aggregate_record(ctx, idx, target, ft, label, model, best)
                if record is not None:
                    records[target][f'{ft}.{label}'] = record
    return records


def emit_aggregates(ctx: RunContext, records: dict) -> None:
    header = aggregate_header(ctx)
    rows = [[f'{t}.{key}'] + [r[k] for k in header[1:]] for t, by_key in records.items() for key, r in by_key.items()]
    write_tsv('run_aggregate_accuracy_mqc.tsv', header, rows)
    plot = {}
    for target, by_key in records.items():
        for record in by_key.values():
            if record['feature_type'] == 'lnc_RNA':
                plot.setdefault(target, {})[record['model']] = round(record['f1'], 2)
    write_json(
        'run_aggregate_accuracy_mqc.json',
        {
            'id': 'run_aggregate_accuracy_plot',
            'section_name': 'allModels vs best single donor — lncRNA F1 per target',
            'description': (
                'lncRNA transcript-level F1 per target for the best single donor and for '
                'allModels (every source pooled, raw), side by side: does pooling every '
                'donor beat the best single one?'
            ),
            'plot_type': 'bargraph',
            'pconfig': {
                'id': 'run_aggregate_accuracy_bargraph',
                'title': 'FOMO: aggregate model accuracy (lncRNA, transcript level)',
                'ylab': 'F1',
                'stacking': None,
                'cpswitch': False,
            },
            'data': plot,
        },
    )


DONOR_HEADER = ['Sample', 'mean_f1', 'median_f1', 'best_f1', 'n_targets_scored', 'mean_models_projected']


def compute_donors(ctx: RunContext, idx: RecordIndex) -> dict:
    """Per source: lncRNA transcript F1 across the benchmarked targets, and how many models it
    projects per target."""
    records = {}
    for source in ctx.sources:
        scores = [
            cell['f1']
            for (target, src, track), cell in idx.accuracy.items()
            if src == source and track == 'lnc_RNA' and target != source
        ]
        counts = [
            n
            for (target, model, track), n in idx.projected_models.items()
            if model == source and track == 'lnc_RNA' and target != source and n is not None
        ]
        records[source] = {
            'mean_f1': statistics.fmean(scores) if scores else None,
            'median_f1': statistics.median(scores) if scores else None,
            'best_f1': max(scores) if scores else None,
            'n_targets_scored': len(scores),
            'mean_models_projected': statistics.fmean(counts) if counts else None,
        }
    return records


# ── entry point ─────────────────────────────────────────────────────────────────
def bool_arg(value: str) -> bool:
    return str(value).lower() in ('true', '1', 'yes')


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--records', required=True, type=Path, help='directory holding the run\'s *.data.json records')
    p.add_argument('--roles', required=True, type=Path, help='species_roles.csv (species,role,has_gff3)')
    p.add_argument('--include-mrna', type=bool_arg, default=False)
    p.add_argument('--include-decoy', type=bool_arg, default=False)
    p.add_argument('--psauron-min', type=float, default=0.5)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not args.records.is_dir():
        fomo_io.die(f"--records '{args.records}' is not a directory")
    feature_types = ('lnc_RNA',) + (('mRNA',) if args.include_mrna else ())
    ctx = RunContext(read_roles(args.roles), feature_types, args.include_decoy, args.psauron_min)
    idx = RecordIndex(args.records)

    metrics, overview = compute_overview(ctx)
    write_tsv('run_overview_mqc.tsv', ['Metric', 'Value'], metrics)
    species = compute_species(ctx, idx)
    write_tsv(
        'run_species_mqc.tsv', SPECIES_HEADER, [[s] + [r[k] for k in SPECIES_HEADER[1:]] for s, r in species.items()]
    )
    funnel, funnel_plot = compute_source_funnel(ctx, idx)
    emit_source_funnel(funnel, funnel_plot)
    targets = compute_targets(ctx, idx)
    emit_targets(ctx, targets)
    heatmap = compute_heatmap(ctx, idx)
    emit_heatmap(ctx, heatmap)
    aggregates = compute_aggregates(ctx, idx)
    emit_aggregates(ctx, aggregates)
    donors = compute_donors(ctx, idx)
    write_tsv(
        'run_donor_ranking_mqc.tsv', DONOR_HEADER, [[s] + [r[k] for k in DONOR_HEADER[1:]] for s, r in donors.items()]
    )

    accuracy_f1 = {
        source: dict(zip(ctx.targets, row, strict=True)) for source, row in zip(ctx.sources, heatmap['f1'], strict=True)
    }
    summary = {
        'overview': overview,
        'sources': ctx.sources,
        'targets': ctx.targets,
        'species': species,
        'source_funnel': funnel,
        'targets_before_after': {t: dict(v) for t, v in targets.items()},
        'accuracy_f1': accuracy_f1,
        'aggregates': {t: dict(v) for t, v in aggregates.items()},
        'donors': donors,
    }
    write_json('run_summary.json', summary)
    return 0


if __name__ == '__main__':
    fomo_io.run(main)
