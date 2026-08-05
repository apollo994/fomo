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
    ch_projected_real   // [ meta(id=source, target_id, feature_type, decoy:false), gff3 ] × 2N  (per-source REAL only)
    ch_all_stats        // [ meta, *.gffcompare.stats ] — per-source stats for ranking (selector filters to real lncRNA)
    ch_target_refs      // [ meta(feature_type), gff3 ] × 2  (filtered target references, from BENCHMARKING)

    main:

    // ── Rank sources by lncRNA transcript F1, keep the top N ──────────────────
    SELECT_TOP_SOURCES(
        ch_all_stats
            .map { _meta, stats -> stats }
            .collect()
            .map { stats -> tuple([id: 'top_sources'], stats) }
    )

    // CSV → channel of selected source ids
    ch_selected = SELECT_TOP_SOURCES.out.csv
        .map { _meta, csv -> csv }
        .splitCsv(header: true)
        .map { row -> tuple(row.source, true) }

    // Keep only the selected sources' projections; group by gene type; combine.
    ch_combine_in = ch_projected_real
        .map { meta, gff -> tuple(meta.id, meta, gff) }
        .combine(ch_selected, by: 0)                      // inner join on source id
        .map { _src, meta, gff, _sel ->
            def gtype = meta.feature_type                 // lnc_RNA | mRNA (no decoys here)
            tuple([id: 'top3', target_id: meta.target_id,
                   feature_type: gtype, decoy: false, gtype: gtype], gff)
        }
        .groupTuple()
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
    ch_top_paired = ch_top_gff3
        .map { meta, gff -> tuple(meta.feature_type, meta, gff) }
        .combine(
            ch_target_refs.map { meta, gff -> tuple(meta.feature_type, gff) },
            by: 0
        )

    ch_top_split = ch_top_paired.multiMap { _ft, q_meta, q_gff, r_gff ->
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

    emit:
    gff3      = ch_top_gff3              // [ meta(id:top3), gff3 ] × 2
    stats     = GFFCOMPARE_TOP.out.stats // [ meta, *.stats ] × 2 (feeds the accuracy scatter in REPORTING)
    mqc_files = ch_mqc_files             // [ meta, path    ] × 4 (2 gffcompare stats + 2 GFF stats tables)
}
