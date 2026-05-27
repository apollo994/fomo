include { PREPROCESSING      } from '../subworkflows/local/preprocessing'
include { PROJECTION         } from '../subworkflows/local/projection'
include { BENCHMARKING       } from '../subworkflows/local/benchmarking'
include { REPORTING          } from '../subworkflows/local/reporting'
include { samplesheetToList  } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { paramsSummaryLog   } from 'plugin/nf-schema'

workflow FOMO {
    // Validate params against nextflow_schema.json (aborts on bad/unknown
    // params) and log a summary of the non-default params in use.
    validateParameters()
    log.info paramsSummaryLog(workflow)

    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .branch { meta, fasta, gff3 ->
            target: meta.role == 'target'
            source: meta.role == 'source'
        }
        .set { ch_input }

    PREPROCESSING(ch_input.source)

    PROJECTION(
        ch_input.target,
        PREPROCESSING.out.spliced_fasta,
        PREPROCESSING.out.decoy_spliced_fasta
    )

    BENCHMARKING(ch_input.target, PROJECTION.out.gff3, PROJECTION.out.target_fasta)

    REPORTING(
        PREPROCESSING.out.mqc_files
            .mix(PROJECTION.out.mqc_files)
            .mix(BENCHMARKING.out.mqc_files)
    )
}
