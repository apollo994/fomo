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

    // ── Roles ────────────────────────────────────────────────────────────────
    // 'both' means donor AND target — the all-vs-all case, where one run replaces
    // N runs with N hand-written samplesheets. A row therefore belongs to two
    // channels, which rules out .branch (it routes each item to exactly ONE
    // output): two .filter's over the same channel, which Nextflow forks.
    //
    // Read the sheet into a plain List first so the three invariants the JSON
    // schema cannot express are checked once, at launch, with a clear message —
    // rather than surfacing later as an empty channel and a pipeline that
    // silently does nothing.
    def rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")

    def dups = rows.collect { meta, _fa, _gff -> meta.id }
                   .countBy { it }.findAll { _id, n -> n > 1 }.keySet()
    if (dups) {
        error("Duplicate species in samplesheet: ${dups.join(', ')}. " +
              "The species id keys every join in the pipeline and must be unique — " +
              "a species that is both donor and target belongs on ONE row with role 'both'.")
    }
    if (!rows.any { meta, _fa, _gff -> meta.role in ['target', 'both'] }) {
        error("Samplesheet has no row with role 'target' or 'both' — nothing to annotate.")
    }
    if (!rows.any { meta, _fa, _gff -> meta.role in ['source', 'both'] }) {
        error("Samplesheet has no row with role 'source' or 'both' — no annotation to transfer.")
    }

    ch_input   = Channel.fromList(rows)
    ch_sources = ch_input.filter { meta, _fasta, _gff3 -> meta.role in ['source', 'both'] }
    ch_targets = ch_input.filter { meta, _fasta, _gff3 -> meta.role in ['target', 'both'] }

    // Source-side work is entirely target-independent (decoys relocate into the
    // SOURCE's own intergenic space), so this runs once per species — S tasks, not
    // S·T. Do not fan it out per target.
    PREPROCESSING(ch_sources, feature_types)

    PROJECTION(
        ch_targets,
        PREPROCESSING.out.spliced_fasta,
        PREPROCESSING.out.decoy_spliced_fasta
    )

    // Only targets carrying a gff3 are benchmarked; BENCHMARKING filters internally
    // and its projection ⋈ reference join is an inner join on [target_id,
    // feature_type], so projections onto an un-annotated target are dropped there.
    BENCHMARKING(ch_targets, PROJECTION.out.gff3, feature_types)

    // Per-source REAL projections only (exclude the all-source 'combined'
    // consensus and decoys) — the pool from which the top-3 consensus is built.
    // Self-pairs are still in here; CONSENSUS_TOP drops them (see below).
    ch_projected_real = PROJECTION.out.gff3
        .filter { meta, _gff -> !meta.decoy && meta.id != 'combined' }

    CONSENSUS_TOP(
        ch_projected_real,
        BENCHMARKING.out.stats,
        BENCHMARKING.out.target_refs
    )

    // One MultiQC report per target. Source-side stats are target-agnostic and get
    // broadcast into every report, so REPORTING needs the target id list to fan
    // them out with.
    REPORTING(
        PREPROCESSING.out.mqc_files
            .mix(PROJECTION.out.mqc_files)
            .mix(BENCHMARKING.out.mqc_files)
            .mix(CONSENSUS_TOP.out.mqc_files),
        BENCHMARKING.out.stats.mix(CONSENSUS_TOP.out.stats),
        ch_targets.map { meta, _fasta, _gff3 -> meta.id }
    )
}
