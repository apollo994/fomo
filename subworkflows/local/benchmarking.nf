include { GFFCOMPARE_BATCH                         } from '../../modules/local/gffcompare_batch'
include { GFF_STATS         as GFF_STATS_TARGET    } from '../../modules/local/gff_stats'
include { pairBy                                   } from './fomo_utils'

workflow BENCHMARKING {
    take:
    ch_target          // [ meta(role:'target'|'both'), fasta[.gz], gff3[.gz]|[] ] × T
    ch_projected_gff3  // [ meta(target_id, track), List<gff3>, manifest ] × F·D·T:
                       // every source's projected GFF3 + the allModels aggregate
    ch_annotation      // [ meta(id, role, feature_type, decoy:false), filtered.gff3 ] × F·S:
                       // each source's filtered annotation (PREPROCESSING.out.filtered_gff3)

    main:

    // The reference is the species' own filtered annotation (spliced, longest isoform,
    // lncRNA TD2-noncoding), so only donors can be benchmarked: the 'both' rows. It is
    // TD2-filtered like the projected models, so coding lncRNA that the pipeline can never
    // report do not count as misses.
    ch_target_refs = ch_annotation.filter { meta, _gff -> meta.role in ['target', 'both'] }

    // A pure target's samplesheet gff3 gets its own raw stats; a 'both' species' annotation is
    // measured on the source side (the same file). target_id routes them to that target's report.
    ch_target_gff = ch_target
        .filter { meta, _fa, gff3 -> meta.role == 'target' && gff3 }
        .map    { meta, _fa, gff3 -> tuple(meta + [kind: 'raw', target_id: meta.id], gff3, []) }

    // Pair each target's batch with the reference of the same target and feature type. The
    // key must include target_id: on feature_type alone a projection would also be compared
    // with every other target's reference, under the same output name. The inner join also
    // drops targets without a reference.
    GFFCOMPARE_BATCH(
        pairBy(
            ch_projected_gff3,
            ch_target_refs.map { meta, gff3 -> tuple([meta.id, meta.feature_type], gff3) },
            { meta -> [meta.target_id, meta.feature_type] }
        )
    )

    GFF_STATS_TARGET(ch_target_gff)

    emit:
    mqc  = GFFCOMPARE_BATCH.out.stats                                // [ meta(target_id), List<*.gffcompare.stats> ]: MultiQC's native gffcompare module
    data = GFFCOMPARE_BATCH.out.data.mix(GFF_STATS_TARGET.out.data)  // [ meta, *.data.json | List ]
}
