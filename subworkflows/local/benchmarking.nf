include { FILTER_TRANSCRIPT as FILTER_TARGET       } from '../../modules/local/filter_transcript'
include { GFFCOMPARE_BATCH                         } from '../../modules/local/gffcompare_batch'
include { GFF_STATS         as GFF_STATS_TARGET    } from '../../modules/local/gff_stats'
include { GFF_STATS_TO_MQC  as GFF_STATS_TARGET_TO_MQC } from '../../modules/local/gff_stats_to_mqc'

workflow BENCHMARKING {
    take:
    ch_target          // [ meta(role:'target'|'both'), fasta[.gz], gff3[.gz]|[] ] × T
    ch_projected_gff3  // [ meta(target_id, feature_type, decoy), List<gff3> ] × F·D·T — per
                       // target: every source's projected GFF3 + the 2 aggregate models
                       // (PROJECTION.out.gff3, see plans/18_per_target_batching.md)
    feature_types      // plain List<String> — MUST be the same list PREPROCESSING got,
                       // hence derived once in workflows/fomo.nf. The projection ⋈
                       // reference pairing below is an INNER join on
                       // [target_id, feature_type], so a list that disagrees with the
                       // source side silently drops projections rather than failing.

    main:

    // Benchmarking needs a reference annotation, so only targets that brought a
    // gff3 take part. A pure 'target' row may omit it (the assembly is annotated
    // but not scored): nf-schema hands a missing path through as [], which is
    // falsy. Its projections are dropped by the inner join below, so no further
    // guard is needed — and no gffcompare, target GFF stats, or top-N consensus
    // are produced for it.

    // The samplesheet gff3 is handed straight to both consumers, still gzipped —
    // there is deliberately NO gunzip task here. FILTER_TARGET decompresses inline
    // (modules/local/filter_transcript.nf:20-24) and gff-feature-stats reads
    // .gff3.gz natively (modules/local/gff_stats.nf:5-8), exactly as
    // PREPROCESSING already relies on for the raw source GFF3. Both derive their
    // output names from ext.prefix, so nothing published depends on the file's
    // name. The ~1 s of decompression now happens once per consumer instead of
    // once per target, which buys back T_g SLURM jobs whose median lifetime was
    // 73 s for 1.1 s of work (see plans/17_scheduling_efficiency.md).
    ch_target_gff = ch_target
        .filter { _meta, _fa, gff3 -> gff3 }
        .map    { meta,  _fa, gff3 -> tuple(meta, gff3) }

    // Filter target annotation by enabled feature_type — one GFF3 per class.
    ch_target_filter_in = ch_target_gff
        .flatMap { meta, gff3 ->
            feature_types.collect { ft -> tuple(meta + [feature_type: ft, decoy: false], gff3, ft) }
        }

    FILTER_TARGET(ch_target_filter_in)

    // Pair each target's whole batch (every source + the 2 aggregates, one list) with
    // the reference of ITS OWN target and matching feature_type. The key MUST include
    // target_id: keyed on feature_type alone every projection would also be compared
    // against every OTHER target's annotation, and since ext.prefix (conf/modules.config)
    // is built from the QUERY meta, all T of those tasks would emit the same filename
    // into the same published directory — wrong numbers under a plausible-looking name.
    //
    // Being an inner join, this is also what drops projections onto targets with
    // no gff3. Unlike the pre-plans/18 per-source combine, this is now a 1:1 join — one
    // batch per (target, feature_type, decoy) against one reference per (target,
    // feature_type) — not a fan-out, since decoy doesn't appear in FILTER_TARGET's key
    // and both decoy states of a target join the same non-decoy reference, same as before.
    ch_paired = ch_projected_gff3
        .map { meta, gffs -> tuple(meta.target_id, meta.feature_type, meta, gffs) }
        .combine(
            FILTER_TARGET.out.gff3.map { meta, gff3 -> tuple(meta.id, meta.feature_type, gff3) },
            by: [0, 1]
        )

    ch_split = ch_paired.multiMap { _tid, _ftype, q_meta, q_gffs, r_gff ->
        query:     tuple(q_meta, q_gffs)
        reference: tuple([id: "${q_meta.target_id}.${q_meta.feature_type}"], r_gff)
    }

    GFFCOMPARE_BATCH(
        ch_split.query,
        ch_split.reference
    )

    // Re-key each stats file back to its own source (or aggregate pseudo-source) — the
    // exact same string-slice re-keying projection.nf uses for the GFF3 split (the split
    // template and this arithmetic must agree; see CLAUDE.md's "merge/split contract").
    // CONSENSUS_TOP:SELECT_TOP_SOURCES needs this: its ranking filter excludes self and
    // the aggregate pseudo-sources by `meta.id` BEFORE grouping per target
    // (consensus_top.nf), so GFFCOMPARE_BATCH's one-shared-meta batch output must be
    // exploded back into individual (meta, statsFile) pairs first. This is a THIRD place
    // that must agree with consensus_top.nf's channel filter and fomo_stats.py's
    // AGGREGATE_IDS (CLAUDE.md's meta-map contract already names the first two).
    ch_stats_persource = GFFCOMPARE_BATCH.out.stats
        .transpose()
        .map { meta, stats ->
            def pre = "${meta.target_id}.from_"
            def suf = ".${meta.feature_type}${meta.decoy ? '.decoy' : ''}.gffcompare.stats"
            tuple([id:           stats.name[pre.size()..-(suf.size() + 1)],
                   feature_type: meta.feature_type,
                   decoy:        meta.decoy,
                   target_id:    meta.target_id], stats)
        }

    // ── Target statistics & MultiQC adapters ─────────────────────────────────
    // target_id is stamped here (== meta.id, the target species) so REPORTING can
    // route these rows into that target's report and no other. Everything else
    // that is target-scoped already carries the key from PROJECTION.
    ch_target_stats_gff = ch_target_gff
        .map  { m, g -> tuple(m + [kind: 'raw',      target_id: m.id], g) }
        .mix(FILTER_TARGET.out.gff3.map { m, g -> tuple(m + [kind: 'filtered', target_id: m.id], g) })

    GFF_STATS_TARGET(ch_target_stats_gff)
    GFF_STATS_TARGET_TO_MQC(GFF_STATS_TARGET.out.json)

    // The custom accuracy scatter is built in REPORTING (over the union of these
    // stats and the top-3 consensus stats) so the consensus shows as its own dot.
    ch_mqc_files = ch_stats_persource
        .mix(GFF_STATS_TARGET_TO_MQC.out.tsv)

    // S = sources, F = enabled feature types, D = 2 with --include_decoy else 1,
    // T_g = targets that supplied a gff3 (everything here scales with T_g, not T).
    // GFFCOMPARE_BATCH itself runs F·D·T_g tasks (plans/18), but the re-keyed stats
    // channel below is unchanged in SHAPE from before batching — one item per source
    // (or aggregate) per target, same (F·S·D + 2·F·D)·T_g count of individual files.
    emit:
    stats         = ch_stats_persource          // [ meta, *.stats ] × (F·S·D + 2·F·D)·T_g
    target_refs   = FILTER_TARGET.out.gff3      // [ meta(id=target, feature_type), gff3 ] × F·T_g
    mqc_files     = ch_mqc_files                // [ meta, path    ] (the stats above + (F+1) GFF_STATS pairs per target: transcript + gene table each)
}
