include { PREPROCESSING      } from '../subworkflows/local/preprocessing'
include { PROJECTION         } from '../subworkflows/local/projection'
include { BENCHMARKING       } from '../subworkflows/local/benchmarking'
include { CONSENSUS_TOP      } from '../subworkflows/local/consensus_top'
include { REPORTING          } from '../subworkflows/local/reporting'
include { samplesheetToList  } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { paramsSummaryLog   } from 'plugin/nf-schema'

workflow FOMO {
    // Validate params against nextflow_schema.json (aborts on bad/unknown
    // params) and log a summary of the non-default params in use.
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // Feature types transferred this run. lncRNA is always on; mRNA is opt-in
    // (--include_mrna) because it is a positive control, not a deliverable.
    // PREPROCESSING (source filtering) and BENCHMARKING (target-reference
    // filtering) MUST agree: benchmarking pairs each projection with the target
    // reference of the same feature_type via an inner combine(by: 0), so a
    // mismatch silently DROPS projections instead of failing. Handing both
    // subworkflows the same list makes that agreement structural — do not
    // re-derive it inside either subworkflow.
    //
    // A plain Groovy List, NOT a channel: it is expanded inside a flatMap closure in
    // each subworkflow. Do not wrap it in Channel.value() and .combine() it — combine
    // SPREADS a List-valued channel into the tuple, so the closure receives the bare
    // String 'lnc_RNA' and String.collect{} then iterates its 7 characters, silently
    // fanning out 7× per source (l/n/c/_/R/N/A tracks). Verified the hard way.
    feature_types = ['lnc_RNA'] + (params.include_mrna ? ['mRNA'] : [])

    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .branch { meta, fasta, gff3 ->
            target: meta.role == 'target'
            source: meta.role == 'source'
        }
        .set { ch_input }

    PREPROCESSING(ch_input.source, feature_types)

    PROJECTION(
        ch_input.target,
        PREPROCESSING.out.spliced_fasta,
        PREPROCESSING.out.decoy_spliced_fasta
    )

    BENCHMARKING(ch_input.target, PROJECTION.out.gff3, feature_types)

    // Per-source REAL projections only (exclude the all-source 'combined'
    // consensus and decoys) — the pool from which the top-3 consensus is built.
    ch_projected_real = PROJECTION.out.gff3
        .filter { meta, _gff -> !meta.decoy && meta.id != 'combined' }

    CONSENSUS_TOP(
        ch_projected_real,
        BENCHMARKING.out.stats,
        BENCHMARKING.out.target_refs
    )

    REPORTING(
        PREPROCESSING.out.mqc_files
            .mix(PROJECTION.out.mqc_files)
            .mix(BENCHMARKING.out.mqc_files)
            .mix(CONSENSUS_TOP.out.mqc_files),
        BENCHMARKING.out.stats.mix(CONSENSUS_TOP.out.stats)
    )
}
