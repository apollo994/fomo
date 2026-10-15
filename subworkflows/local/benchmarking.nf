include { GFFCOMPARE_BATCH                         } from '../../modules/local/gffcompare_batch'
include { GFF_STATS         as GFF_STATS_TARGET    } from '../../modules/local/gff_stats'
include { GFF_STATS_TO_MQC  as GFF_STATS_TARGET_TO_MQC } from '../../modules/local/gff_stats_to_mqc'

workflow BENCHMARKING {
    take:
    ch_target          // [ meta(role:'target'|'both'), fasta[.gz], gff3[.gz]|[] ] × T
    ch_projected_gff3  // [ meta(target_id, feature_type, decoy), List<gff3> ] × F·D·T:
                       // every source's projected GFF3 + the allModels aggregate
    ch_annotation      // [ meta(id, role, feature_type, decoy:false), filtered.gff3 ] × F·S:
                       // each source's filtered annotation (PREPROCESSING.out.filtered_gff3)

    main:

    // The reference is the species' own filtered annotation (spliced, longest isoform,
    // lncRNA TD2-noncoding), so only donors can be benchmarked: the 'both' rows. It is
    // TD2-filtered like the projected models, so coding lncRNA that the pipeline can never
    // report do not count as misses.
    ch_target_refs = ch_annotation.filter { meta, _gff -> meta.role in ['target', 'both'] }

    // Every target with a samplesheet gff3 gets its raw annotation stats.
    ch_target_gff = ch_target
        .filter { _meta, _fa, gff3 -> gff3 }
        .map    { meta,  _fa, gff3 -> tuple(meta, gff3) }

    // Pair each target's batch with the reference of the same target and feature type. The
    // key must include target_id: on feature_type alone a projection would also be compared
    // with every other target's reference, under the same output name. The inner join also
    // drops targets without a reference.
    ch_paired = ch_projected_gff3
        .map { meta, gffs -> tuple(meta.target_id, meta.feature_type, meta, gffs) }
        .combine(
            ch_target_refs.map { meta, gff3 -> tuple(meta.id, meta.feature_type, gff3) },
            by: [0, 1]
        )

    GFFCOMPARE_BATCH(ch_paired.map { _tid, _ftype, meta, gffs, ref -> tuple(meta, gffs, ref) })

    // One item per stats file; consumers read the target from meta and the source from the
    // file name.
    ch_stats_persource = GFFCOMPARE_BATCH.out.stats.transpose()

    // target_id (== the species) routes these rows to that target's report only.
    ch_target_stats_gff = ch_target_gff
        .map  { m, g -> tuple(m + [kind: 'raw',      target_id: m.id], g) }
        .mix(ch_target_refs.map { m, g -> tuple(m + [kind: 'filtered', target_id: m.id], g) })

    GFF_STATS_TARGET(ch_target_stats_gff)
    GFF_STATS_TARGET_TO_MQC(GFF_STATS_TARGET.out.json)

    ch_mqc_files = ch_stats_persource
        .mix(GFF_STATS_TARGET_TO_MQC.out.tsv)

    emit:
    stats     = ch_stats_persource   // [ meta(target_id, feature_type, decoy), *.gffcompare.stats ]
    mqc_files = ch_mqc_files         // [ meta, path ]
}
