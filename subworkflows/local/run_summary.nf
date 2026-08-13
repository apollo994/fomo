include { RUN_SUMMARY_TABLES        } from '../../modules/local/run_summary_tables'
include { MULTIQC as MULTIQC_RUN_SUMMARY } from '../../modules/nf-core/multiqc/main'

// ONE run-level report, alongside REPORTING's one-per-target ones. It answers the
// questions no per-target report can: how many species in which role, how much
// annotation went in and how much survived filtering, how much each target gained,
// and which donors are worth using.
//
// Aggregate-ONLY by design. The obvious alternative — hand the whole union to a
// second MULTIQC — produces a report BIGGER than the per-target ones (at S=T=20 the
// union is ~1000 files / ~900 table rows) and the per-target accuracy scatters all
// carry the same hardcoded section id (`gffcompare_accuracy_custom`), so they would
// overwrite each other. One task digests the union into a dozen run-level tables
// instead.
//
// Deliberately a TOP-LEVEL subworkflow, not nested and not folded into REPORTING:
// every selector in conf/modules.config is a fully-qualified process path and
// `withName` matching is a regex *find*, so `FOMO:REPORTING:MULTIQC` would also
// match a `MULTIQC_RUN_SUMMARY` alias declared inside REPORTING and hand it the
// per-target ext.prefix. A separate subworkflow makes that miss structural.
workflow RUN_SUMMARY {
    take:
    ch_all_mqc      // Channel<tuple(meta, path|List<path>)> — the pipeline-wide union
    ch_top_sources  // Channel<tuple(meta, *.top_sources.csv)> × T_g
    ch_roles        // Channel<path> — species_roles.csv (single file)

    main:

    // flatMap, NOT map: GFF_STATS_TO_MQC emits a LIST of two paths per item (the
    // transcript table and the gene-category table) — the same trap REPORTING
    // documents. Everything is then collected into ONE task; the script selects
    // files by suffix out of the staged directory, so no per-suffix channel
    // splitting is needed here.
    ch_files = ch_all_mqc
        .flatMap { _meta, p -> p instanceof List ? p : [p] }
        .collect()

    // ifEmpty([]) is load-bearing: with no annotated target there is no ranking and
    // no top_sources.csv, and .collect() on an empty channel emits NOTHING — the
    // process would wait forever on an input that never arrives.
    ch_top = ch_top_sources
        .map { _meta, csv -> csv }
        .collect()
        .ifEmpty([])

    RUN_SUMMARY_TABLES(ch_files, ch_top, ch_roles)

    MULTIQC_RUN_SUMMARY(
        RUN_SUMMARY_TABLES.out.mqc.map { files ->
            tuple(
                // No target_id, on purpose: this report belongs to the run, not to a
                // target. Its publishDir rule is a static path for the same reason —
                // a target-scoped closure here would publish to targets/null/.
                [id: 'run_summary'],
                files,
                [
                    file(params.multiqc_run_summary_main_config,     checkIfExists: true),
                    file(params.multiqc_run_summary_sections_config, checkIfExists: true)
                ],
                [],
                [],
                []
            )
        }
    )

    emit:
    report = MULTIQC_RUN_SUMMARY.out.report   // [ meta(id:'run_summary'), fomo_run_summary.html ]
    tables = RUN_SUMMARY_TABLES.out.tables    // run_summary.json — every number, machine-readable
}
