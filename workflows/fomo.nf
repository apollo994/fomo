include { PREPROCESSING      } from '../subworkflows/local/preprocessing'
include { PROJECTION         } from '../subworkflows/local/projection'
include { BENCHMARKING       } from '../subworkflows/local/benchmarking'
include { CURATION           } from '../subworkflows/local/curation'
include { REPORTING          } from '../subworkflows/local/reporting'
include { RUN_SUMMARY        } from '../subworkflows/local/run_summary'
include { softwareVersionsYaml } from '../subworkflows/local/fomo_utils'
include { samplesheetToList  } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { paramsSummaryLog   } from 'plugin/nf-schema'

workflow FOMO {
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // Feature types transferred this run: lncRNA always, mRNA (a positive control) on request.
    // Derived here only: BENCHMARKING's references come from PREPROCESSING, already split by
    // this list, so the projection ⋈ reference join cannot disagree with it. A plain List,
    // not a channel (see PREPROCESSING).
    feature_types = ['lnc_RNA'] + (params.include_mrna ? ['mRNA'] : [])

    // A 'both' row is a donor and a target, so it goes to both channels: two filters, not
    // branch. The schema rejects duplicate species; the checks below cover what it cannot.
    def rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")

    if (!rows.any { meta, _fa, _gff -> meta.role in ['target', 'both'] }) {
        error("Samplesheet has no row with role 'target' or 'both' — nothing to annotate.")
    }
    if (!rows.any { meta, _fa, _gff -> meta.role in ['source', 'both'] }) {
        error("Samplesheet has no row with role 'source' or 'both' — no annotation to transfer.")
    }

    ch_input   = channel.fromList(rows)
    ch_sources = ch_input.filter { meta, _fasta, _gff3 -> meta.role in ['source', 'both'] }
    ch_targets = ch_input.filter { meta, _fasta, _gff3 -> meta.role in ['target', 'both'] }

    // A pure target's gff3 is its curation reference only: benchmarking needs the filtered
    // annotation that only donors have.
    def ref_only = rows.findAll { meta, _fa, gff3 -> meta.role == 'target' && gff3 }.collect { meta, _fa, _gff -> meta.id }
    if (ref_only) {
        log.info "Targets with a reference annotation used for curation only (not benchmarked — " +
                 "use role 'both' to benchmark, which also makes them donors): ${ref_only.join(', ')}"
    }

    PREPROCESSING(ch_sources, feature_types)

    PROJECTION(
        ch_targets,
        PREPROCESSING.out.spliced_fasta,
        PREPROCESSING.out.decoy_spliced_fasta
    )

    BENCHMARKING(ch_targets, PROJECTION.out.gff3, PREPROCESSING.out.filtered_gff3)

    CURATION(PROJECTION.out.allmodels_raw, ch_targets)

    // The MultiQC inputs of the per-target reports, and the data records of the run summary.
    ch_all_mqc = PREPROCESSING.out.mqc_files
        .mix(PROJECTION.out.mqc_files)
        .mix(BENCHMARKING.out.mqc_files)
        .mix(CURATION.out.mqc_files)

    ch_target_ids = ch_targets.map { meta, _fasta, _gff3 -> meta.id }

    // Every module puts its tool versions on the `versions` topic. The YAML goes to the MultiQC
    // tasks only: a consumer that also emits versions would wait on itself.
    ch_versions = softwareVersionsYaml(channel.topic('versions'))

    REPORTING(
        ch_all_mqc,
        BENCHMARKING.out.stats,
        ch_target_ids,
        ch_versions
    )

    // Species roles are the one fact no stats file carries; the sorted CSV keeps the task cacheable.
    ch_roles = channel.fromList(rows)
        .map { meta, _fasta, gff3 -> "${meta.id},${meta.role},${gff3 ? 'yes' : 'no'}" }
        .collectFile(name: 'species_roles.csv', newLine: true, sort: true,
                     seed: 'species,role,has_gff3')

    ch_all_data = PREPROCESSING.out.data
        .mix(PROJECTION.out.data)
        .mix(BENCHMARKING.out.data)
        .mix(CURATION.out.data)

    RUN_SUMMARY(ch_all_data, ch_roles, ch_versions)
}
