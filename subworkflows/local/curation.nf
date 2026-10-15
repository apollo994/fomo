include { CURATE_MODELS } from '../../modules/local/curate_models'

// Curated annotation, per (target, feature type, decoy): intron chains shared exactly by
// >= params.curate_min_species species, minus models with an exon on a same-strand reference
// exon, one transcript per supported chain, merged into the target's reference.
// Rationale: plans/21, plans/25.
workflow CURATION {
    take:
    ch_allmodels_raw    // [ meta(id:'allModels', target_id, feature_type, decoy), allModels.<gtype>.raw.gff3 ] × F·D·T
    ch_targets          // [ meta(id = target, role), fasta, gff3 | [] ] × T

    main:

    // The reference is the target's samplesheet GFF3 (every biotype), not the filtered
    // annotation: overlaps with coding genes, tRNAs and pseudogenes must count. A target
    // without a gff3 is curated without a reference (no filter, no merged file). Keyed on
    // target_id alone: one reference serves every track of a target.
    ch_refs = ch_targets.map { meta, _fasta, gff3 -> tuple(meta.id, gff3 ?: []) }

    ch_curate_in = ch_allmodels_raw
        .map { meta, gff -> tuple(meta.target_id, meta, gff) }
        .combine(ch_refs, by: 0)
        .map { _tid, meta, gff, ref -> tuple(meta + [id: 'curated'], gff, ref) }

    CURATE_MODELS(ch_curate_in)

    emit:
    mqc_files = CURATE_MODELS.out.mqc   // [ meta(target_id), [ 3 × *_mqc.tsv ] ]: that target's report only
}
