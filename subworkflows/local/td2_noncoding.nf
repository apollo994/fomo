include { TD2_PREDICT       } from '../../modules/local/td2_predict'
include { TD2_CODING_FILTER } from '../../modules/local/td2_coding_filter'

// Run TD2 on a lncRNA FASTA and split coding from non-coding. A transcript is
// "coding" if TD2.Predict finds a complete ORF with psauron_score above the
// configured threshold (see TD2_CODING_FILTER ext.args). Reused for both the
// input lncRNA (preprocessing) and the per-source projected lncRNA (projection).
workflow TD2_NONCODING {
    take:
    ch_fasta   // [ meta, fasta ]  — lncRNA transcript FASTA

    main:
    TD2_PREDICT(ch_fasta)

    // Pair each FASTA with its own .pep by join on meta (same meta object flows
    // through TD2_PREDICT, so a positional join is safe and order-independent).
    ch_filter_in = ch_fasta.join(TD2_PREDICT.out.pep)

    TD2_CODING_FILTER(ch_filter_in)

    emit:
    kept_fasta = TD2_CODING_FILTER.out.kept_fasta   // [ meta, *.noncoding.fasta ]
    coding_ids = TD2_CODING_FILTER.out.coding_ids   // [ meta, *.coding_ids.txt  ]
    mqc        = TD2_CODING_FILTER.out.mqc          // [ meta, *_td2_coding_mqc.tsv ]
}
