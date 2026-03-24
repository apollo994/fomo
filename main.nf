#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO_BUSCO    } from './modules/nf-core/busco/busco/main'
include { BUSCO_DOWNLOAD } from './modules/local/busco_download/main'

workflow RUN_BUSCO {
    main:
    Channel
        .fromPath(params.input, checkIfExists: true)
        .map { fasta -> [ [id: fasta.baseName], fasta ] }
        .set { ch_fasta }

    BUSCO_DOWNLOAD(params.lineage)

    BUSCO_BUSCO(
        ch_fasta,
        params.busco_mode,
        params.lineage,
        BUSCO_DOWNLOAD.out.lineage_path,
        [],
        false
    )

    emit:
    batch_summary = BUSCO_BUSCO.out.batch_summary
    busco_dir     = BUSCO_BUSCO.out.busco_dir
}

workflow {
    RUN_BUSCO()
}
