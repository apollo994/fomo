include { SELECT_TOP_SOURCES                       } from '../../modules/local/select_top_sources'
include { GFF_BY_SOURCE as SUBSET_GFF_BY_SOURCE    } from '../../modules/local/gff_by_source'
include { GFFCOMPARE as GFFCOMPARE_COMBINE_TOP     } from '../../modules/nf-core/gffcompare/main'
include { GFFREAD    as COMBINED_TOP_GTF_TO_GFF    } from '../../modules/nf-core/gffread/main'
include { GFF_STATS         as GFF_STATS_TOP       } from '../../modules/local/gff_stats'
include { GFF_STATS_TO_MQC  as GFF_STATS_TOP_TO_MQC } from '../../modules/local/gff_stats_to_mqc'
include { GFFCOMPARE as GFFCOMPARE_TOP             } from '../../modules/nf-core/gffcompare/main'

// Build a "top-N consensus" annotation: rank sources by lncRNA transcript-level F1
// (from per-source gffcompare stats), pull those sources' models straight out of the
// all-sources annotation, and collapse them with gffcompare. Both the raw subset and
// the collapsed consensus are scored, under the pseudo-source ids 'top3_raw' and
// 'top3_collapsed'.
workflow CONSENSUS_TOP {
    take:
    ch_allmodels_raw    // [ meta(id='allModels', target_id, feature_type, decoy, gtype), gff3 ] × F·D·T
    ch_all_stats        // [ meta, *.gffcompare.stats ] — per-source stats for ranking (filtered below)
    ch_target_refs      // [ meta(id=target, feature_type), gff3 ] × F·T_g  (filtered target references, from BENCHMARKING)

    main:

    // ── Rank sources by lncRNA transcript F1, keep the top N — PER TARGET ─────
    // One SELECT_TOP_SOURCES task per target, not one globally: the ranking is a
    // property of the target (which relatives transfer best onto THIS assembly),
    // and bin/select_top_sources.py keys its score map by source alone, so pooling
    // several targets' stats into one task would have them overwrite each other.
    //
    // Two things are excluded from the pool:
    //
    // * SELF. A species projected onto itself scores near 100% F1 and would always
    //   take the #1 slot, turning the consensus into a restatement of the target's
    //   own annotation. This is the ONLY place self is excluded — it is deliberately
    //   kept in allModels (projection.nf).
    // * The AGGREGATE pseudo-sources. allModels_raw / allModels_collapsed /
    //   top3_raw / top3_collapsed travel through the same `from_<id>` naming closures
    //   as per-source projections, so their stats land in this same channel; leaving
    //   them in would pick the top-N out of models that are themselves built from the
    //   top-N. bin/select_top_sources.py repeats the list as a second layer — keep
    //   the two in sync.
    ch_rank_in = ch_all_stats
        .filter { meta, _stats ->
            !meta.decoy && meta.feature_type == 'lnc_RNA' &&
            !(meta.id in ['allModels_raw', 'allModels_collapsed', 'top3_raw', 'top3_collapsed']) &&
            meta.id != meta.target_id
        }
        .map { meta, stats -> tuple(meta.target_id, stats) }
        .groupTuple()
        .map { tid, stats -> tuple([id: tid, target_id: tid], stats) }

    SELECT_TOP_SOURCES(ch_rank_in)

    // ── Extract the selected sources' models from allModels ───────────────────
    // A subset of the same file the per-source models were split from, so top3.raw
    // is exactly "allModels restricted to the top N sources" — no re-derivation, no
    // possibility of drift between the two. Decoys have no ranking and are excluded.
    //
    // combine(by: 0) on target_id is an inner join, which is also what drops targets
    // with no gff3: they produce no ranking stats, hence no CSV, hence nothing here.
    // Joining on the target — not broadcasting the CSV — is what stops target A's
    // top-N being applied to target B's models.
    ch_subset_in = ch_allmodels_raw
        .filter { meta, _gff -> !meta.decoy }
        .map { meta, gff -> tuple(meta.target_id, meta, gff) }
        .combine(
            SELECT_TOP_SOURCES.out.csv.map { meta, csv -> tuple(meta.target_id, csv) },
            by: 0
        )
        .map { _tid, meta, gff, csv -> tuple(meta + [id: 'top3_raw'], gff, csv) }

    SUBSET_GFF_BY_SOURCE(ch_subset_in)

    ch_top_raw = SUBSET_GFF_BY_SOURCE.out.gff3.transpose()

    // ── Collapse the subset into the consensus ────────────────────────────────
    GFFCOMPARE_COMBINE_TOP(
        ch_top_raw,
        [[:], [], []],   // no reference sequence  (-s)
        [[:], []]        // no reference annotation (-r) → pure combine mode
    )

    COMBINED_TOP_GTF_TO_GFF(
        GFFCOMPARE_COMBINE_TOP.out.combined_gtf.map { meta, gtf -> tuple(meta + [id: 'top3_collapsed'], gtf) },
        []   // no genome FASTA needed for a pure GTF→GFF3 conversion
    )

    // Both models are scored: the raw subset answers "how good is the union of the
    // best N sources?", the collapsed one "…and what does merging it cost?".
    ch_top_scored = ch_top_raw.mix(COMBINED_TOP_GTF_TO_GFF.out.gffread_gff)

    // ── GFF statistics on both consensus models ───────────────────────────────
    GFF_STATS_TOP(
        ch_top_scored.map { meta, gff -> tuple(meta + [kind: 'projected'], gff) }
    )
    GFF_STATS_TOP_TO_MQC(GFF_STATS_TOP.out.json)

    // ── gffcompare each against the matching target reference ─────────────────
    // Keyed on [target_id, feature_type] for the same reason as the equivalent
    // join in benchmarking.nf — on feature_type alone each consensus would also be
    // scored against every other target's annotation, under a filename built from
    // the query meta and therefore identical across all T results.
    ch_top_paired = ch_top_scored
        .map { meta, gff -> tuple(meta.target_id, meta.feature_type, meta, gff) }
        .combine(
            ch_target_refs.map { meta, gff -> tuple(meta.id, meta.feature_type, gff) },
            by: [0, 1]
        )

    ch_top_split = ch_top_paired.multiMap { _tid, _ft, q_meta, q_gff, r_gff ->
        query:     tuple(q_meta, q_gff)
        empty_ref: tuple([id: q_meta.target_id], [], [])
        reference: tuple([id: "${q_meta.target_id}.${q_meta.feature_type}"], r_gff)
    }

    GFFCOMPARE_TOP(
        ch_top_split.query,
        ch_top_split.empty_ref,
        ch_top_split.reference
    )

    ch_mqc_files = GFF_STATS_TOP_TO_MQC.out.tsv
        .mix(GFFCOMPARE_TOP.out.stats)

    // F = enabled feature types, T_g = targets that supplied a gff3. A target with
    // no annotation contributes no ranking stats, so it produces nothing here —
    // no guard needed, the empty join simply never matches.
    emit:
    gff3        = ch_top_scored            // [ meta(id:top3_raw|top3_collapsed), gff3 ] × 2·F·T_g
    stats       = GFFCOMPARE_TOP.out.stats // [ meta, *.stats ] × 2·F·T_g (feeds the accuracy scatter in REPORTING)
    mqc_files   = ch_mqc_files             // [ meta, path    ] × 4·F·T_g (gffcompare stats + GFF stats tables)
    // The selection itself, for RUN_SUMMARY. Not a MultiQC file (one bare `source`
    // column), and the F1 ranking behind it is only printed to stderr, so this CSV
    // is the only record of which donors a target actually chose. EMPTY when no
    // target is annotated — consumers must .ifEmpty([]) before collecting.
    top_sources = SELECT_TOP_SOURCES.out.csv // [ meta(id=target), <target>.top_sources.csv ] × T_g
}
