include { CURATE_MODELS } from '../../modules/local/curate_models'

// Curated annotation (plans/21): per (target, feature type, decoy) — every enabled track,
// lncRNA, mRNA and their decoys — keep only intron chains shared exactly by
// >= params.curate_min_species species, drop models overlapping any reference exon (same
// strand), group supported chains sharing a splice junction into genes, report one
// representative per gene, and merge the result into the target's reference annotation.
//
// A TOP-LEVEL subworkflow (called from workflows/fomo.nf, not nested in PROJECTION):
// every `withName` selector in conf/modules.config is a fully-qualified path
// (FOMO:CURATION:CURATE_MODELS) — see CLAUDE.md, multi-target rule 1.
workflow CURATION {
    take:
    ch_allmodels_raw    // [ meta(id:'allModels', target_id, feature_type, decoy, gtype), allModels.<gtype>.raw.gff3 ] × F·D·T
    ch_targets          // [ meta(id = target, role), fasta, gff3 | [] ] × T

    main:

    // The reference is the target's SAMPLESHEET GFF3 — the full annotation, every biotype —
    // not PREPROCESSING's filtered ★ annotation (lncRNA-only, longest isoform, TD2-filtered):
    // "ignore models overlapping the reference" must see coding genes, tRNAs, pseudogenes, …
    // A pure target has no gff3 ([] after samplesheetToList) and is curated without a
    // reference — no reference filter, no merged file.
    //
    // combine(by: 0) on target_id: every target is exactly one ch_targets row, so each
    // allModels track meets exactly one reference slot. Keyed on target_id ALONE on
    // purpose — the reference is per target, shared by every feature type and decoy track.
    ch_refs = ch_targets.map { meta, _fasta, gff3 -> tuple(meta.id, gff3 ?: []) }

    ch_curate_in = ch_allmodels_raw
        .map { meta, gff -> tuple(meta.target_id, meta, gff) }
        .combine(ch_refs, by: 0)
        .map { _tid, meta, gff, ref -> tuple(meta + [id: 'curated'], gff, ref) }

    CURATE_MODELS(ch_curate_in)

    emit:
    curated   = CURATE_MODELS.out.curated   // [ meta(id:'curated'), <target>.curated.<gtype>.gff3.gz ] × F·D·T
    merged    = CURATE_MODELS.out.merged    // [ meta, <target>.curated.<gtype>.merged.gff3.gz ] × (tracks of annotated targets)
    report    = CURATE_MODELS.out.report    // [ meta, [ .curation.tsv, .curation.json ] ]
    mqc_files = CURATE_MODELS.out.mqc       // [ meta(target_id), [ 3 × *_mqc.tsv ] ] → that target's report only
}
