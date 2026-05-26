include { MULTIQC } from '../../modules/nf-core/multiqc/main'

workflow REPORTING {
    take:
    ch_mqc_files   // Channel<tuple(meta, path)> — *_mqc.tsv / *.stats files

    main:

    ch_files = ch_mqc_files.map { _meta, p -> p }.collect()

    ch_mqc_input = ch_files.map { files ->
        tuple(
            [id: 'multiqc'],
            files,
            [
                file(params.multiqc_main_config,     checkIfExists: true),
                file(params.multiqc_sections_config, checkIfExists: true)
            ],
            [],
            [],
            []
        )
    }

    MULTIQC(ch_mqc_input)

    emit:
    report = MULTIQC.out.report
    data   = MULTIQC.out.data
}
