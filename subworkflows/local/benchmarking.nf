include { GUNZIP as GUNZIP_TARGET_GFF        } from '../../modules/nf-core/gunzip/main'
include { FILTER_TRANSCRIPT as FILTER_TARGET } from '../../modules/local/filter_transcript'
include { GFFCOMPARE                         } from '../../modules/nf-core/gffcompare/main'

workflow BENCHMARKING {
    take:
    ch_target          // [ meta(role:'target'), fasta.gz, gff3.gz ] × 1
    ch_projected_gff3  // [ meta(target_id, id, feature_type, decoy), gff3 ] × 16

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

    emit:
    stats = GFFCOMPARE.out.stats    // [ meta, *.stats ] × 16
}
