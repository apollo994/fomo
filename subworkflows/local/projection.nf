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
    ch_target            // [ meta(role:'target'), fasta[.gz], gff3[.gz] ] × 1
    ch_spliced_fasta     // [ meta(decoy:false, feature_type), fasta ] × 2N
    ch_decoy_spliced     // [ meta(decoy:true,  feature_type), fasta ] × 2N

    main:

    // Decompress the target FASTA once
    GUNZIP_TARGET(
        ch_target.map { meta, fa, _gff3 -> tuple(meta, fa) }
    )

    // Build .mmi index once (ext.args = '-x splice' set in conf/modules.config)
    MINIMAP2_INDEX(GUNZIP_TARGET.out.gunzip)

    // Broadcast index as value channel so it's reused across every alignment task
    ch_reference = MINIMAP2_INDEX.out.index
        .map { meta, mmi -> tuple([id: meta.id], mmi) }
        .first()

    // Reads channel: every source spliced FASTA (mRNA / lnc_RNA / decoy).
    // Attach target_id to meta so output filenames encode both target and source.
    ch_reads = ch_spliced_fasta
        .mix(ch_decoy_spliced)
        .combine(ch_reference)
        .map { meta, reads, ref_meta, _mmi ->
            tuple(meta + [target_id: ref_meta.id], reads)
        }

    MINIMAP2_ALIGN(
        ch_reads,
        ch_reference,
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

    // Single target genome FASTA, broadcast to every gffread task.
    ch_target_fa = GUNZIP_TARGET.out.gunzip.map { _m, fa -> fa }.first()

    EXTRACT_PROJECTED_LNC(ch_proj.lnc, ch_target_fa)

    TD2_NONCODING_PROJ(EXTRACT_PROJECTED_LNC.out.gffread_fasta)

    // Drop coding lncRNA from each per-source projected GFF3 (join gff + ids on
    // the shared meta, preserved unchanged through gffread + TD2).
    ch_proj_lnc_filter_in = ch_proj.lnc.join(TD2_NONCODING_PROJ.out.coding_ids)
    FILTER_PROJ_LNC_GFF(ch_proj_lnc_filter_in)

    // Filtered real lncRNA ⊕ untouched (real mRNA + all decoys).
    ch_projected_persource = FILTER_PROJ_LNC_GFF.out.gff3.mix(ch_proj.rest)

    // ── Combine projected models per gene type ────────────────────────────────
    // Group every projected GFF3 by gene type — feature_type plus the decoy
    // flag → {lnc_RNA, mRNA, lnc_RNA_decoy, mRNA_decoy} — and run gffcompare in
    // combine mode (no reference annotation) to build one consensus
    // <prefix>.combined.gtf per gene type across all source species.
    ch_combine_in = ch_projected_persource
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

    emit:
    bam          = MINIMAP2_ALIGN.out.bam            // [ meta, *.bam          ]
    index        = MINIMAP2_ALIGN.out.index          // [ meta, *.bam.bai      ]
    gff3         = ch_projected                      // [ meta, *.gff3         ] × 20 (16 per-source + 4 combined)
    combined_gff = COMBINED_GTF_TO_GFF.out.gffread_gff // [ meta, *.gff3        ] × 4 (consensus per gene type)
    mqc_files    = ch_mqc_files                      // [ meta, *_mqc.tsv      ] (one table per projected model + samtools + TD2)
}
