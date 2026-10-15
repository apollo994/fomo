include { TD2_PREDICT       } from '../../modules/local/td2_predict'
include { TD2_CODING_FILTER } from '../../modules/local/td2_coding_filter'

// Split a lncRNA FASTA into coding and non-coding: a transcript is coding when TD2.Predict
// finds a complete ORF with psauron_score >= params.td2_psauron_min. Used on the source lncRNA
// (PREPROCESSING) and on each target's pooled projected lncRNA (PROJECTION).
workflow TD2_NONCODING {
    take:
    ch_fasta   // [ meta(stage: 'input' | 'projected'), fasta ]: lncRNA transcript FASTA

    main:
    TD2_PREDICT(ch_fasta)

    ch_filter_in = ch_fasta.join(TD2_PREDICT.out.pep)

    TD2_CODING_FILTER(ch_filter_in)

    emit:
    kept_fasta = TD2_CODING_FILTER.out.kept_fasta   // [ meta, *.noncoding.fasta ]
    coding_ids = TD2_CODING_FILTER.out.coding_ids   // [ meta, *.coding_ids.txt  ]
    mqc        = TD2_CODING_FILTER.out.mqc          // [ meta, *_td2_coding_mqc.tsv ]
    data       = TD2_CODING_FILTER.out.data         // [ meta, *.td2_coding.data.json ]
}
