include { PREPROCESSING      } from '../subworkflows/local/preprocessing'
include { PROJECTION         } from '../subworkflows/local/projection'
include { BENCHMARKING       } from '../subworkflows/local/benchmarking'
include { CURATION           } from '../subworkflows/local/curation'
include { REPORTING          } from '../subworkflows/local/reporting'
include { RUN_SUMMARY        } from '../subworkflows/local/run_summary'
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
    // Consumed by PREPROCESSING only. BENCHMARKING's references are PREPROCESSING's
    // own filtered annotations (plans/20_longest_isoform.md), already split by this
    // list, so the projection ⋈ reference inner join cannot disagree with the source
    // side on feature_type — do not re-derive the list inside any subworkflow.
    //
    // A plain Groovy List, NOT a channel: it is expanded inside a flatMap closure in
    // PREPROCESSING. Do not wrap it in Channel.value() and .combine() it — combine
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

    // A pure 'target' may carry a gff3 (plans/22): it is that target's CURATION reference
    // (filter overlapping models, tag the nearest reference gene, merged output) and its raw
    // GFF stats are reported — but it is NOT benchmarked, since the gffcompare reference is the
    // filtered annotation that only donors have (plans/20). Say so at launch, so a missing
    // gffcompare section for such a target is not a surprise.
    def ref_only = rows.findAll { meta, _fa, gff3 -> meta.role == 'target' && gff3 }.collect { meta, _fa, _gff -> meta.id }
    if (ref_only) {
        log.info "Targets with a reference annotation used for curation only (not benchmarked — " +
                 "use role 'both' to benchmark, which also makes them donors): ${ref_only.join(', ')}"
    }

    // Source-side work is entirely target-independent (decoys relocate into the
    // SOURCE's own intergenic space), so this runs once per species — S tasks, not
    // S·T. Do not fan it out per target.
    PREPROCESSING(ch_sources, feature_types)

    PROJECTION(
        ch_targets,
        PREPROCESSING.out.spliced_fasta,
        PREPROCESSING.out.decoy_spliced_fasta
    )

    // Only targets carrying a gff3 are benchmarked — the 'both' rows, since the schema
    // rejects a gff3 on a pure 'target'. The reference is the species' ONE filtered
    // annotation from PREPROCESSING (the file its own projection was extracted from);
    // the projection ⋈ reference join is an inner join on [target_id, feature_type],
    // so projections onto an un-annotated target are dropped there.
    BENCHMARKING(ch_targets, PROJECTION.out.gff3, PREPROCESSING.out.filtered_gff3)


    // Curated annotation (plans/21): every allModels track — lncRNA, mRNA, decoys — keeps
    // only intron chains shared exactly by >= params.curate_min_species species, minus any
    // model overlapping the target's reference, and is merged into that reference. The
    // reference is the samplesheet GFF3 (ch_targets), not the filtered ★ annotation.
    // `--curate false` is applied as `ext.when` on CURATE_MODELS (conf/modules.config), not
    // an `if` here: the process then always exists, so its config selectors never warn.
    CURATION(PROJECTION.out.allmodels_raw, ch_targets)

    // The pipeline-wide union of MultiQC inputs, hoisted because it feeds BOTH
    // reporting paths. A channel read by two consumers is forked by Nextflow, so
    // this is a naming change only — neither consumer sees a shortened stream.
    ch_all_mqc = PREPROCESSING.out.mqc_files
        .mix(PROJECTION.out.mqc_files)
        .mix(BENCHMARKING.out.mqc_files)
        .mix(CURATION.out.mqc_files)

    ch_target_ids = ch_targets.map { meta, _fasta, _gff3 -> meta.id }

    // One MultiQC report per target. Source-side stats are target-agnostic and get
    // broadcast into every report, so REPORTING needs the target id list to fan
    // them out with.
    REPORTING(
        ch_all_mqc,
        BENCHMARKING.out.stats,
        ch_target_ids
    )

    // ...and ONE report for the run as a whole, which is the only place the
    // cross-target picture exists: species per role, how much annotation survived
    // filtering, what each target gained, which donors are worth using.
    //
    // Species roles are the one fact no stat file carries, so they are passed as a
    // small CSV built from the samplesheet rows read above. `seed` is written first,
    // `sort` makes the body deterministic (and therefore the task cacheable).
    ch_roles = Channel.fromList(rows)
        .map { meta, _fasta, gff3 -> "${meta.id},${meta.role},${gff3 ? 'yes' : 'no'}" }
        .collectFile(name: 'species_roles.csv', newLine: true, sort: true,
                     seed: 'species,role,has_gff3')

    RUN_SUMMARY(ch_all_mqc, ch_roles)
}
