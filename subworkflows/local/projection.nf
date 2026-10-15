include { MERGE_FASTA as MERGE_SOURCE_FASTA        } from '../../modules/local/merge_fasta'
include { MINIMAP2_INDEX                           } from '../../modules/local/minimap2_index'
include { MINIMAP2_ALIGN                           } from '../../modules/local/minimap2_align'
include { SAMTOOLS_STATS                           } from '../../modules/local/samtools_stats'
include { SAMTOOLS_TO_MQC                          } from '../../modules/local/samtools_to_mqc'
include { BAM_TO_GFF                               } from '../../modules/local/bam_to_gff'
include { GFF_STATS_PROJECTED_BATCH                } from '../../modules/local/gff_stats_projected_batch'
include { EXTRACT_TRANSCRIPTS as EXTRACT_PROJECTED_LNC } from '../../modules/local/extract_transcripts'
include { TD2_NONCODING     as TD2_NONCODING_PROJ  } from './td2_noncoding'
include { FILTER_GFF_BY_ID  as FILTER_ALLMODELS    } from '../../modules/local/filter_gff_by_id'
include { GFF_BY_SOURCE     as SPLIT_GFF_BY_SOURCE } from '../../modules/local/gff_by_source'

workflow PROJECTION {
    take:
    ch_target            // [ meta(role:'target'|'both'), fasta[.gz], gff3[.gz]|[] ] × T
    ch_spliced_fasta     // [ meta(decoy:false, feature_type), fasta ] × F·S
    ch_decoy_spliced     // [ meta(decoy:true,  feature_type), fasta ] × F·S

    main:

    MINIMAP2_INDEX(ch_target.map { meta, fa, _gff3 -> tuple(meta, fa) })

    // [ target_id, is_source, mmi ]: a 'both' target's own records are dropped at alignment.
    ch_ref = MINIMAP2_INDEX.out.index.map { meta, mmi -> tuple(meta.id, meta.role == 'both', mmi) }

    // One all-sources FASTA per gene type (F·D tasks, target-independent). Every record carries
    // its species (`>tid|type|species`), which the split below recovers.
    ch_merged = ch_spliced_fasta
        .mix(ch_decoy_spliced)
        .map { meta, fa ->
            def gtype = meta.feature_type + (meta.decoy ? '_decoy' : '')
            tuple([id: gtype, feature_type: meta.feature_type, decoy: meta.decoy, gtype: gtype], fa)
        }
        .groupTuple()
        .map { meta, fas -> tuple(meta, fas.sort { fa -> fa.name }) }   // deterministic concat order

    MERGE_SOURCE_FASTA(ch_merged)

    // One alignment per (target, gene type): F·D·T, each against its target's index.
    // exclude_self drops the target's own records from the query, so a species is never
    // projected onto itself.
    MINIMAP2_ALIGN(
        MERGE_SOURCE_FASTA.out.fasta
            .combine(ch_ref)
            .map { meta, reads, tid, is_source, mmi ->
                tuple(meta + [id: 'allModels', target_id: tid, exclude_self: is_source], reads, mmi)
            }
    )

    SAMTOOLS_STATS(MINIMAP2_ALIGN.out.bam)

    // A curated column set for the report; the native samtools section keeps its mapped-%
    // chart. One row per (target, gene type): the merged BAM has no per-source breakdown.
    SAMTOOLS_TO_MQC(SAMTOOLS_STATS.out.stats)

    BAM_TO_GFF(MINIMAP2_ALIGN.out.bam.map { meta, bam, _csi -> tuple(meta, bam) })

    // TD2 runs once per target on the pooled lncRNA models; the other tracks pass untouched.
    BAM_TO_GFF.out.gff3
        .branch { meta, _g ->
            lnc:  !meta.decoy && meta.feature_type == 'lnc_RNA'
            rest: true
        }
        .set { ch_proj }

    // Transcripts are extracted from the target's own genome (joined on target_id).
    ch_target_fa = ch_target.map { meta, fa, _gff3 -> tuple(meta.id, fa) }

    EXTRACT_PROJECTED_LNC(
        ch_proj.lnc
            .map { meta, gff -> tuple(meta.target_id, meta, gff) }
            .combine(ch_target_fa, by: 0)
            .map { _tid, meta, gff, fa -> tuple(meta, gff, fa) }
    )

    TD2_NONCODING_PROJ(EXTRACT_PROJECTED_LNC.out.fasta.map { meta, fa -> tuple(meta + [stage: 'projected'], fa) })

    // Every track goes through FILTER_ALLMODELS, so one process writes every
    // `allModels.<gtype>.raw.gff3`; only lncRNA has a drop list (one per target).
    ch_filter_in = ch_proj.lnc
        .map { meta, gff -> tuple(meta.target_id, meta, gff) }
        .join(TD2_NONCODING_PROJ.out.coding_ids.map { meta, ids -> tuple(meta.target_id, ids) })
        .map { _tid, meta, gff, ids -> tuple(meta, gff, ids, []) }
        .mix(ch_proj.rest.map { meta, gff -> tuple(meta, gff, [], []) })

    FILTER_ALLMODELS(ch_filter_in)

    ch_allmodels_raw = FILTER_ALLMODELS.out.gff3   // [ meta(id:'allModels'), gff3 ] × F·D·T

    // Back to one GFF3 per source, by the `source=` attribute.
    SPLIT_GFF_BY_SOURCE(ch_allmodels_raw)

    // One list per (target, gene type): the aggregate plus every source's file. A target onto
    // which nothing projected has no split output; `remainder: true` keeps it, so it is
    // reported with zero models instead of disappearing.
    ch_projected_batch = ch_allmodels_raw
        .join(SPLIT_GFF_BY_SOURCE.out.gff3, remainder: true)
        .map { meta, raw_gff, split_gffs -> tuple(meta, ([raw_gff] + [split_gffs].flatten()).findAll { gff -> gff }) }

    GFF_STATS_PROJECTED_BATCH(
        ch_projected_batch.map { meta, gffs -> tuple(meta + [kind: 'projected'], gffs) }
    )

    ch_mqc_files = GFF_STATS_PROJECTED_BATCH.out.tsv
        .mix(SAMTOOLS_STATS.out.stats)
        .mix(SAMTOOLS_TO_MQC.out.tsv)
        .mix(TD2_NONCODING_PROJ.out.mqc)

    emit:
    allmodels_raw = ch_allmodels_raw     // [ meta(id:'allModels'), *.gff3 ] × F·D·T
    gff3          = ch_projected_batch   // [ meta, List<gff3> ] × F·D·T: allModels + every source
    mqc_files     = ch_mqc_files         // [ meta, path | List<path> ]
}
