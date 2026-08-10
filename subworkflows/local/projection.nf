include { MAYBE_GUNZIP as GUNZIP_TARGET            } from '../../modules/local/maybe_gunzip'
include { MINIMAP2_INDEX                           } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                           } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_STATS                           } from '../../modules/nf-core/samtools/stats/main.nf'
include { SAMTOOLS_TO_MQC                          } from '../../modules/local/samtools_to_mqc'
include { BAM_TO_GFF                               } from '../../modules/local/bam_to_gff'
include { GFF_STATS         as GFF_STATS_PROJECTED } from '../../modules/local/gff_stats'
include { GFF_STATS_TO_MQC  as GFF_STATS_PROJECTED_TO_MQC } from '../../modules/local/gff_stats_to_mqc'
include { GFFCOMPARE        as GFFCOMPARE_COMBINE  } from '../../modules/nf-core/gffcompare/main'
include { GFFREAD           as COMBINED_GTF_TO_GFF } from '../../modules/nf-core/gffread/main'
include { GFFREAD           as EXTRACT_PROJECTED_LNC } from '../../modules/nf-core/gffread/main'
include { TD2_NONCODING     as TD2_NONCODING_PROJ  } from './td2_noncoding'
include { FILTER_GFF_BY_ID  as FILTER_PROJ_LNC_GFF } from '../../modules/local/filter_gff_by_id'

workflow PROJECTION {
    take:
    ch_target            // [ meta(role:'target'|'both'), fasta[.gz], gff3[.gz]|[] ] × T
    ch_spliced_fasta     // [ meta(decoy:false, feature_type), fasta ] × F·S
    ch_decoy_spliced     // [ meta(decoy:true,  feature_type), fasta ] × F·S

    main:

    // Decompress each target FASTA once (T tasks)
    GUNZIP_TARGET(
        ch_target.map { meta, fa, _gff3 -> tuple(meta, fa) }
    )

    // Build one .mmi index per target (ext.args = '-x splice' set in conf/modules.config)
    MINIMAP2_INDEX(GUNZIP_TARGET.out.gunzip)

    // Per-target reference bundle, keyed by target id: [ target_id, mmi, fasta ].
    // MINIMAP2_INDEX and GUNZIP_TARGET both preserve the input meta, so join() on
    // the full meta pairs each index with its own genome.
    ch_ref = MINIMAP2_INDEX.out.index
        .join(GUNZIP_TARGET.out.gunzip)
        .map { meta, mmi, fa -> tuple(meta.id, mmi, fa) }

    // The all-vs-all fan-out: every source track (mRNA / lnc_RNA / decoy) × every
    // target. combine() without `by` is a full cartesian product and buffers the
    // right-hand side, so each of the T references is reused across all F·S·D
    // tracks. Self-pairs (source == target) are KEPT — they are the Sn/Pr ceiling
    // control — but excluded from the consensus further down.
    //
    // multiMap, not a value channel: with T > 1 the reference must be paired with
    // its own reads task-by-task. A broadcast/.first() reference would silently
    // align every source onto whichever target's index materialised first.
    ch_align = ch_spliced_fasta
        .mix(ch_decoy_spliced)
        .combine(ch_ref)
        .multiMap { meta, reads, tid, mmi, _fa ->
            reads:     tuple(meta + [target_id: tid], reads)
            reference: tuple([id: tid], mmi)
        }

    MINIMAP2_ALIGN(
        ch_align.reads,
        ch_align.reference,
        true,    // bam_format
        'bai',   // bam_index_extension
        false,   // cigar_paf_format
        false    // cigar_bam
    )

    // Join BAM + BAI into a single 3-element tuple [meta, bam, bai]
    ch_bam_bai = MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index)

    SAMTOOLS_STATS(
        ch_bam_bai,
        [[:], [], []]       // no reference FASTA needed for basic alignment stats
    )

    // Parse `*.stats` into a curated 7-column MultiQC custom-content TSV that
    // surfaces as a section under Samtools (MultiQC's native samtools-stats
    // module hard-codes which columns appear in its Alignment Stats violin
    // plot and can't be customised via config).
    SAMTOOLS_TO_MQC(SAMTOOLS_STATS.out.stats)

    BAM_TO_GFF(MINIMAP2_ALIGN.out.bam)

    // ── TD2 coding-potential filter on projected lncRNA (real, per-source) ────
    // Extract each per-source projected lncRNA's transcript FASTA from the
    // target genome, run TD2, and drop coding-potential models from the GFF3
    // BEFORE the per-gene-type combine — so combined + top-3 consensus inherit
    // the filtering. Real mRNA and all decoys pass through untouched.
    BAM_TO_GFF.out.gff3
        .branch { meta, _g ->
            lnc:  !meta.decoy && meta.feature_type == 'lnc_RNA'
            rest: true
        }
        .set { ch_proj }

    // Each projection's transcript FASTA must be extracted from ITS OWN target
    // genome — hence the join on target_id rather than a broadcast FASTA. Getting
    // this wrong is silent: gffread would happily read the coordinates of target B
    // out of target A's sequence and hand TD2 nonsense to score.
    ch_target_fa = GUNZIP_TARGET.out.gunzip.map { meta, fa -> tuple(meta.id, fa) }

    ch_lnc_extract = ch_proj.lnc
        .map { meta, gff -> tuple(meta.target_id, meta, gff) }
        .combine(ch_target_fa, by: 0)
        .multiMap { _tid, meta, gff, fa ->
            gff:   tuple(meta, gff)
            fasta: fa
        }

    EXTRACT_PROJECTED_LNC(ch_lnc_extract.gff, ch_lnc_extract.fasta)

    TD2_NONCODING_PROJ(EXTRACT_PROJECTED_LNC.out.gffread_fasta)

    // Drop coding lncRNA from each per-source projected GFF3 (join gff + ids on
    // the shared meta, preserved unchanged through gffread + TD2).
    ch_proj_lnc_filter_in = ch_proj.lnc.join(TD2_NONCODING_PROJ.out.coding_ids)
    FILTER_PROJ_LNC_GFF(ch_proj_lnc_filter_in)

    // Filtered real lncRNA ⊕ untouched (real mRNA + all decoys).
    ch_projected_persource = FILTER_PROJ_LNC_GFF.out.gff3.mix(ch_proj.rest)

    // ── Combine projected models per (target, gene type) ──────────────────────
    // Group every projected GFF3 by target and gene type — feature_type plus the
    // decoy flag → {lnc_RNA, mRNA, lnc_RNA_decoy, mRNA_decoy} — and run gffcompare
    // in combine mode (no reference annotation) to build one consensus
    // <target>.<gtype>.combined.gtf per gene type across all source species.
    //
    // Self-pairs are EXCLUDED here: a species' own annotation projected onto itself
    // is a near-perfect copy, so leaving it in would make the consensus mostly a
    // restatement of the target's existing annotation instead of evidence
    // transferred from relatives. It stays available as a per-source model and as
    // the ceiling dot in the accuracy scatter. Consequence: a target whose only
    // source is itself yields no combined model.
    ch_combine_in = ch_projected_persource
        .filter { meta, _gff -> meta.id != meta.target_id }
        .map { meta, gff ->
            def gtype = meta.feature_type + (meta.decoy ? '_decoy' : '')
            tuple([id: "${meta.target_id}.${gtype}", target_id: meta.target_id,
                   feature_type: meta.feature_type, decoy: meta.decoy, gtype: gtype], gff)
        }
        .groupTuple()
        .map { meta, gffs -> tuple(meta, gffs.sort { it.name }) }   // deterministic input order

    GFFCOMPARE_COMBINE(
        ch_combine_in,
        [[:], [], []],   // no reference sequence  (-s)
        [[:], []]        // no reference annotation (-r) → pure combine mode
    )

    // Convert each combined GTF → GFF3 so it matches the rest of the pipeline,
    // and reshape meta to the projection convention (id = 'combined') so the
    // consensus models flow through the same GFF stats / benchmarking naming closures
    // as per-source projections (yielding "from_combined" rows).
    COMBINED_GTF_TO_GFF(
        GFFCOMPARE_COMBINE.out.combined_gtf.map { meta, gtf -> tuple(meta + [id: 'combined'], gtf) },
        []   // no genome FASTA needed for a pure GTF→GFF3 conversion
    )

    // Union of per-source projections (lncRNA coding-filtered) + per-gene-type
    // consensus. All downstream stats/benchmarking treat the combined models as
    // an extra "from_combined" source.
    ch_projected = ch_projected_persource.mix(COMBINED_GTF_TO_GFF.out.gffread_gff)

    // ── Statistics on projected models ───────────────────────────────────────
    GFF_STATS_PROJECTED(
        ch_projected.map { meta, gff -> tuple(meta + [kind: 'projected'], gff) }
    )
    GFF_STATS_PROJECTED_TO_MQC(GFF_STATS_PROJECTED.out.json)

    // Keep SAMTOOLS_STATS.out.stats in the mix so MultiQC's native samtools
    // module still runs (it renders the "Percent mapped" bar chart). Our
    // curated 7-column TSV is rendered as a custom section nested under
    // Samtools (see assets/multiqc/sections.yml: parent_id: samtools).
    ch_mqc_files = GFF_STATS_PROJECTED_TO_MQC.out.tsv
        .mix(SAMTOOLS_STATS.out.stats)
        .mix(SAMTOOLS_TO_MQC.out.tsv)
        .mix(TD2_NONCODING_PROJ.out.mqc)

    // T = targets, S = sources, F = enabled feature types, D = 2 with --include_decoy else 1.
    emit:
    bam          = MINIMAP2_ALIGN.out.bam            // [ meta, *.bam     ] × F·S·T·D
    index        = MINIMAP2_ALIGN.out.index          // [ meta, *.bam.bai ] × F·S·T·D
    gff3         = ch_projected                      // [ meta, *.gff3    ] × F·S·T·D + F·T·D (per-source + combined)
    combined_gff = COMBINED_GTF_TO_GFF.out.gffread_gff // [ meta, *.gff3  ] × F·T·D (consensus per target × gene type)
    mqc_files    = ch_mqc_files                      // [ meta, *_mqc.tsv ] (one table per projected model + samtools + TD2)
}
