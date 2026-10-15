include { MULTIQC       } from '../../modules/local/multiqc'
include { routeByTarget } from './fomo_utils'

// One MultiQC report per target. MULTIQC renders its custom content from the records; native
// tool files feed MultiQC's own modules. A record or file whose meta has a target_id belongs to
// that target's report; the others (source side) go to every report.
workflow REPORTING {
    take:
    ch_records      // [ meta, *.data.json | List ]: every data record of the run
    ch_native       // [ meta, path | List ]: native MultiQC inputs (samtools and gffcompare .stats)
    ch_target_ids   // String: one id per target
    ch_versions     // path: the software versions YAML, shared by every report

    main:

    ch_native_all = ch_native.mix(ch_versions.map { yml -> tuple([:], yml) })

    ch_mqc_input = routeByTarget(ch_records, ch_target_ids)
        .join(routeByTarget(ch_native_all, ch_target_ids), remainder: true)
        .map { tid, records, tool_files ->
            tuple(
                [id: tid, target_id: tid],
                records ?: [],
                tool_files ?: [],
                [
                    file(params.multiqc_main_config,     checkIfExists: true),
                    file(params.multiqc_sections_config, checkIfExists: true)
                ]
            )
        }

    MULTIQC(ch_mqc_input)
}
