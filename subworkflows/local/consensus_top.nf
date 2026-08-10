include { SELECT_TOP_SOURCES                       } from '../../modules/local/select_top_sources'
include { GFFCOMPARE as GFFCOMPARE_COMBINE_TOP     } from '../../modules/nf-core/gffcompare/main'
include { GFFREAD    as COMBINED_TOP_GTF_TO_GFF    } from '../../modules/nf-core/gffread/main'
include { GFF_STATS         as GFF_STATS_TOP       } from '../../modules/local/gff_stats'
include { GFF_STATS_TO_MQC  as GFF_STATS_TOP_TO_MQC } from '../../modules/local/gff_stats_to_mqc'
include { GFFCOMPARE as GFFCOMPARE_TOP             } from '../../modules/nf-core/gffcompare/main'

// Build a "top-N consensus" annotation: rank sources by lncRNA transcript-level
// F1 (from per-source gffcompare stats), then gffcompare-combine the lncRNA AND
// mRNA projections of those same top sources. The consensus is plugged into the
// same downstream as every other model (GFF stats, gffcompare-vs-reference, MultiQC)
// under the pseudo-source id 'top3'.
workflow CONSENSUS_TOP {
    take:
    ch_projected_real   // [ meta(id=source, target_id, feature_type, decoy:false), gff3 ] × F·S·T  (per-source REAL only)
    ch_all_stats        // [ meta, *.gffcompare.stats ] — per-source stats for ranking (filtered below)
    ch_target_refs      // [ meta(id=target, feature_type), gff3 ] × F·T_g  (filtered target references, from BENCHMARKING)

    main:

    // ── Rank sources by lncRNA transcript F1, keep the top N — PER TARGET ─────
    // One SELECT_TOP_SOURCES task per target, not one globally: the ranking is a
    // property of the target (which relatives transfer best onto THIS assembly),
    // and bin/select_top_sources.py keys its score map by source alone, so pooling
    // several targets' stats into one task would have them overwrite each other.
    //
    // Self is dropped from the ranking pool here. A species projected onto itself
    // scores near 100% F1 and would always take the #1 slot, turning the consensus
    // into a restatement of the target's own annotation. Filtering in the channel
    // (rather than in the script) keeps the script target-agnostic — it already
    // excludes decoys and the 'combined'/'top3' pseudo-sources by filename.
    ch_rank_in = ch_all_stats
        .filter { meta, _stats ->
            !meta.decoy && meta.feature_type == 'lnc_RNA' &&
            !(meta.id in ['combined', 'top3']) && meta.id != meta.target_id
        }
        .map { meta, stats -> tuple(meta.target_id, stats) }
        .groupTuple()
        .map { tid, stats -> tuple([id: tid, target_id: tid], stats) }

    SELECT_TOP_SOURCES(ch_rank_in)

    // CSV → [ target_id, source ]. `elem: 1` splits the file while preserving the
    // meta, which is what carries the target this ranking belongs to.
    ch_selected = SELECT_TOP_SOURCES.out.csv
        .splitCsv(header: true, elem: 1)
        .map { meta, row -> tuple(meta.target_id, row.source, true) }

    // Keep only each target's own selected sources; group by (target, gene type);
    // combine. Joining on [target_id, source] — not source alone — is what stops
    // target A's top-3 pulling in projections made onto target B.
    ch_combine_in = ch_projected_real
        .filter { meta, _gff -> meta.id != meta.target_id }   // no self, as above
        .map { meta, gff -> tuple(meta.target_id, meta.id, meta, gff) }
        .combine(ch_selected, by: [0, 1])                    // inner join on (target, source)
        .map { _tid, _src, meta, gff, _sel ->
            def gtype = meta.feature_type                 // lnc_RNA | mRNA (no decoys here)
            tuple([id: 'top3', target_id: meta.target_id,
                   feature_type: gtype, decoy: false, gtype: gtype], gff)
        }
        .groupTuple()                                     // meta carries target_id → per-target groups
        .map { meta, gffs -> tuple(meta, gffs.sort { it.name }) }   // deterministic order

    GFFCOMPARE_COMBINE_TOP(
        ch_combine_in,
        [[:], [], []],   // no reference sequence  (-s)
        [[:], []]        // no reference annotation (-r) → pure combine mode
    )

    // ── Convert consensus GTF → GFF3 (matches the rest of the pipeline) ───────
    COMBINED_TOP_GTF_TO_GFF(
        GFFCOMPARE_COMBINE_TOP.out.combined_gtf,
        []   // no genome FASTA needed for a pure GTF→GFF3 conversion
    )
    ch_top_gff3 = COMBINED_TOP_GTF_TO_GFF.out.gffread_gff   // [ meta(id:top3,...), gff3 ] × 2

    // ── GFF statistics on the consensus ───────────────────────────────────────
    GFF_STATS_TOP(
        ch_top_gff3.map { meta, gff -> tuple(meta + [kind: 'projected'], gff) }
    )
    GFF_STATS_TOP_TO_MQC(GFF_STATS_TOP.out.json)

    // ── gffcompare the consensus against the matching target reference ────────
    // Keyed on [target_id, feature_type] for the same reason as the equivalent
    // join in benchmarking.nf — on feature_type alone each consensus would also be
    // scored against every other target's annotation, under a filename built from
    // the query meta and therefore identical across all T results.
    ch_top_paired = ch_top_gff3
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
    // no guard needed, the empty group simply never forms.
    emit:
    gff3      = ch_top_gff3              // [ meta(id:top3, target_id), gff3 ] × F·T_g
    stats     = GFFCOMPARE_TOP.out.stats // [ meta, *.stats ] × F·T_g (feeds the accuracy scatter in REPORTING)
    mqc_files = ch_mqc_files             // [ meta, path    ] × 2·F·T_g (gffcompare stats + GFF stats tables)
}
