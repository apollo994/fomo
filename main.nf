#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO_BUSCO    } from './modules/nf-core/busco/busco/main'
include { BUSCO_DOWNLOAD } from './modules/local/busco_download/main'

workflow FOMO {
    main:
    Channel
        .fromPath(params.input, checkIfExists: true)
        .map { fasta -> [ [id: fasta.baseName], fasta ] }
        .set { ch_fasta }

    if (params.busco_lineages_path) {
        ch_lineage_path = Channel.value(file(params.busco_lineages_path))
    } else {
        BUSCO_DOWNLOAD(params.lineage)
        ch_lineage_path = BUSCO_DOWNLOAD.out.lineage_path
    }

    BUSCO_BUSCO(
        ch_fasta,
        params.busco_mode,
        params.lineage,
        ch_lineage_path,
        [],
        true
    )

    emit:
    short_summaries_json = BUSCO_BUSCO.out.short_summaries_json
}

workflow {
    FOMO()
}
