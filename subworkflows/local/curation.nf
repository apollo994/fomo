include { CURATE_MODELS } from '../../modules/local/curate_models'
include { pairBy        } from './fomo_utils'

// Curated annotation, per (target, feature type, decoy): intron chains shared exactly by
// >= params.curate_min_species species, minus models with an exon on a same-strand reference
// exon, one transcript per supported chain, merged into the target's reference.
// Rationale: plans/21, plans/25.
workflow CURATION {
    take:
    ch_allmodels_raw    // [ meta(target_id, track), allModels.<track>.raw.gff3 ] × F·D·T
    ch_targets          // [ meta(id = target, role), fasta, gff3 | [] ] × T

    main:

    // The reference is the target's samplesheet GFF3 (every biotype), not the filtered
    // annotation: overlaps with coding genes, tRNAs and pseudogenes must count. A target
    // without a gff3 is curated without a reference (no filter, no merged file). Keyed on
    // target_id alone: one reference serves every track of a target.
    ch_refs = ch_targets.map { meta, _fasta, gff3 -> tuple(meta.id, gff3 ?: []) }

    // `--curate false`: no input, no task.
    ch_models = params.curate ? ch_allmodels_raw : channel.empty()

    CURATE_MODELS(pairBy(ch_models, ch_refs, { meta -> meta.target_id }))

    emit:
    data = CURATE_MODELS.out.data   // [ meta(target_id), *.curation.data.json ]
}
