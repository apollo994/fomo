include { RUN_SUMMARY_TABLES        } from '../../modules/local/run_summary_tables'
include { MULTIQC as MULTIQC_RUN_SUMMARY } from '../../modules/local/multiqc'

// One run-level report: species per role, what survived filtering, what each target gained,
// which donors are worth using. Aggregate-only: one task digests the pipeline-wide union into
// run-level tables, and only those reach its MultiQC. A top-level subworkflow, so the
// per-target `FOMO:REPORTING:MULTIQC$` settings cannot match its MultiQC alias.
workflow RUN_SUMMARY {
    take:
    ch_all_mqc      // [ meta, path | List<path> ]: the pipeline-wide union
    ch_roles        // species_roles.csv
    ch_versions     // path: the software versions YAML

    main:

    // One task; the script selects files by suffix from the staged directory.
    ch_files = ch_all_mqc
        .flatMap { _meta, p -> [p].flatten() }
        .collect()

    RUN_SUMMARY_TABLES(ch_files, ch_roles)

    MULTIQC_RUN_SUMMARY(
        RUN_SUMMARY_TABLES.out.mqc.mix(ch_versions).collect().map { files ->
            tuple(
                // No target_id: the report belongs to the run (static publishDir path).
                [id: 'run_summary'],
                files,
                [
                    file(params.multiqc_run_summary_main_config,     checkIfExists: true),
                    file(params.multiqc_run_summary_sections_config, checkIfExists: true)
                ]
            )
        }
    )
}
