include { GUNZIP as GUNZIP_TARGET_GFF              } from '../../modules/nf-core/gunzip/main'
include { FILTER_TRANSCRIPT as FILTER_TARGET       } from '../../modules/local/filter_transcript'
include { GFFCOMPARE                               } from '../../modules/nf-core/gffcompare/main'
include { AGAT_SPSTATISTICS as AGAT_TARGET         } from '../../modules/local/agat_spstatistics'
include { AGAT_TO_MQC       as AGAT_TARGET_TO_MQC  } from '../../modules/local/agat_to_mqc'

workflow BENCHMARKING {
    take:
    ch_target          // [ meta(role:'target'), fasta.gz, gff3.gz ] × 1
    ch_projected_gff3  // [ meta(target_id, id, feature_type, decoy), gff3 ] × 16
    ch_target_fasta    // [ meta(role:'target'), fasta ] × 1  (gunzipped, from PROJECTION)

    main:

    GUNZIP_TARGET_GFF(
        ch_target.map { meta, _fa, gff3 -> tuple(meta, gff3) }
    )

    // Filter target annotation by feature_type — produces one GFF3 per class.
    ch_target_filter_in = GUNZIP_TARGET_GFF.out.gunzip
        .combine(Channel.of('lnc_RNA', 'mRNA'))
        .map { meta, gff3, ftype ->
            tuple(meta + [feature_type: ftype, decoy: false], gff3, ftype)
        }

    FILTER_TARGET(ch_target_filter_in)

    // Pair each projection with the target-class reference of matching feature_type.
    ch_paired = ch_projected_gff3
        .map { meta, gff3 -> tuple(meta.feature_type, meta, gff3) }
        .combine(
            FILTER_TARGET.out.gff3.map { meta, gff3 -> tuple(meta.feature_type, gff3) },
            by: 0
        )

    ch_split = ch_paired.multiMap { _ftype, q_meta, q_gff, r_gff ->
        query:     tuple(q_meta, q_gff)
        empty_ref: tuple([id: q_meta.target_id], [], [])
        reference: tuple([id: "${q_meta.target_id}.${q_meta.feature_type}"], r_gff)
    }

    GFFCOMPARE(
        ch_split.query,
        ch_split.empty_ref,
        ch_split.reference
    )

    // ── Target statistics & MultiQC adapters ─────────────────────────────────
    ch_target_agat_gff = GUNZIP_TARGET_GFF.out.gunzip
        .map  { m, g -> tuple(m + [kind: 'raw'],      g) }
        .mix(FILTER_TARGET.out.gff3.map { m, g -> tuple(m + [kind: 'filtered'], g) })

    // Pair each target GFF with the (gunzipped) target genome FASTA so AGAT
    // can compute genome-coverage metrics via --gs.
    ch_target_agat_in = ch_target_agat_gff
        .map { meta, gff -> tuple(meta.id, meta, gff) }
        .combine(
            ch_target_fasta.map { meta, fa -> tuple(meta.id, fa) },
            by: 0
        )
        .map { _id, meta, gff, fa -> tuple(meta, gff, fa) }

    AGAT_TARGET(ch_target_agat_in)
    AGAT_TARGET_TO_MQC(AGAT_TARGET.out.stats_yaml)

    ch_mqc_files = GFFCOMPARE.out.stats
        .mix(AGAT_TARGET_TO_MQC.out.tsv)

    emit:
    stats     = GFFCOMPARE.out.stats    // [ meta, *.stats ] × 16
    mqc_files = ch_mqc_files            // [ meta, path    ] × 19 (16 stats + 3 AGAT)
}
