include { GFFCOMPARE_TO_SCATTER } from '../../modules/local/gffcompare_to_scatter'
include { MULTIQC               } from '../../modules/nf-core/multiqc/main'

workflow REPORTING {
    take:
    ch_mqc_files        // Channel<tuple(meta, path)> — *_mqc.tsv / *.stats files
    ch_gffcompare_stats // Channel<tuple(meta, *.gffcompare.stats)> — all benchmark + consensus stats

    main:

    // Build the custom gffcompare accuracy scatter over ALL stats (per-source +
    // all-source combined + top-3 consensus) so every model — including the
    // consensus — gets its own dot. Replaces MultiQC's native accuracy plot
    // (removed via remove_sections in assets/multiqc/main.yml).
    GFFCOMPARE_TO_SCATTER(
        ch_gffcompare_stats
            .map { _meta, stats -> stats }
            .collect()
            .map { stats -> tuple([id: 'gffcompare_accuracy'], stats) }
    )

    ch_files = ch_mqc_files.map { _meta, p -> p }
        .mix(GFFCOMPARE_TO_SCATTER.out.json.map { _meta, p -> p })
        .collect()

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
