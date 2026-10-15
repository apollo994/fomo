include { GFFCOMPARE_TO_SCATTER } from '../../modules/local/gffcompare_to_scatter'
include { MULTIQC               } from '../../modules/nf-core/multiqc/main'

// One MultiQC report per target. Files with meta.target_id belong to that target's report;
// the others (source-side stats) are broadcast into every report.
workflow REPORTING {
    take:
    ch_mqc_files        // [ meta, path | List<path> ]: *_mqc.tsv and native tool stats
    ch_gffcompare_stats // [ meta(target_id), *.gffcompare.stats ]: every source and allModels
    ch_target_ids       // String: one id per target

    main:

    // The accuracy scatter replaces MultiQC's native gffcompare plot (removed in main.yml):
    // one point per model of the target.
    GFFCOMPARE_TO_SCATTER(
        ch_gffcompare_stats
            .map { meta, stats -> tuple(meta.target_id, stats) }
            .groupTuple()
            .map { tid, stats -> tuple([id: tid, target_id: tid], stats) }
    )

    ch_mqc_files
        .branch { meta, _p ->
            scoped: meta.target_id != null
            shared: true
        }
        .set { ch_mqc }

    // One file per item before combine: an adapter may emit a list of files, and combine
    // spreads a list across the tuple.
    ch_by_target = ch_mqc.shared
        .flatMap { _meta, p -> [p].flatten() }
        .combine(ch_target_ids)
        .map { p, tid -> tuple(tid, p) }
        .mix(ch_mqc.scoped.flatMap { meta, p ->
            [p].flatten().collect { f -> tuple(meta.target_id, f) } })
        .mix(GFFCOMPARE_TO_SCATTER.out.json.map { meta, p -> tuple(meta.target_id, p) })
        .groupTuple()

    ch_mqc_input = ch_by_target.map { tid, files ->
        tuple(
            [id: tid, target_id: tid],
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
}
