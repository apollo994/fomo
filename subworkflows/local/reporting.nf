include { GFFCOMPARE_TO_SCATTER } from '../../modules/local/gffcompare_to_scatter'
include { MULTIQC               } from '../../modules/nf-core/multiqc/main'

// One MultiQC report per target. Statistics come in two flavours: target-scoped
// (projections, alignments, gffcompare, the target's own reference stats — all
// carrying meta.target_id) and target-agnostic (the source-side input GFF stats,
// SeqKit, input TD2 — identical whichever target they are read alongside). The
// first go to their own target's report; the second are broadcast into every one.
workflow REPORTING {
    take:
    ch_mqc_files        // Channel<tuple(meta, path)> — *_mqc.tsv / *.stats files
    ch_gffcompare_stats // Channel<tuple(meta, *.gffcompare.stats)> — all benchmark + consensus stats
    ch_target_ids       // Channel<String> — one id per target (role 'target' or 'both')

    main:

    // Build the custom gffcompare accuracy scatter per target, over that target's
    // stats (per-source + all-source combined + top-3 consensus) so every model —
    // including the consensus — gets its own dot. Replaces MultiQC's native
    // accuracy plot (removed via remove_sections in assets/multiqc/main.yml).
    // Within one target's plot the point names are `<source>.<feature_type>` and
    // stay unique, so bin/gffcompare_accuracy_mqc.py needs no change.
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

    // Broadcast the shared files against the target ids, one FILE per item — note
    // the flatMap. Two things force it:
    //   * some adapters emit a LIST of paths per item (GFF_STATS_TO_MQC emits both
    //     a transcript table and a gene-category table), and
    //   * combine SPREADS a List-valued item across the output tuple (the same
    //     trap documented for feature_types in workflows/fomo.nf),
    // so without flattening first, `combine` turns [meta, [a, b]] into [a, b, tid]
    // and the downstream closure is called with three arguments. Equally, do NOT
    // .collect() the shared files into one list and combine THAT in — same spread.
    ch_by_target = ch_mqc.shared
        .flatMap { _meta, p -> p instanceof List ? p : [p] }
        .combine(ch_target_ids)
        .map { p, tid -> tuple(tid, p) }
        .mix(ch_mqc.scoped.flatMap { meta, p ->
            (p instanceof List ? p : [p]).collect { f -> tuple(meta.target_id, f) } })
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

    emit:
    report = MULTIQC.out.report   // [ meta(id=target), *.html ] × T
    data   = MULTIQC.out.data     // [ meta(id=target), *_data ] × T
}
