include { MAYBE_GUNZIP as GUNZIP_TARGET            } from '../../modules/local/maybe_gunzip'
include { MINIMAP2_INDEX                           } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                           } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_STATS                           } from '../../modules/nf-core/samtools/stats/main.nf'
include { SAMTOOLS_TO_MQC                          } from '../../modules/local/samtools_to_mqc'
include { BAM_TO_GFF                               } from '../../modules/local/bam_to_gff'
include { AGAT_SPSTATISTICS as AGAT_PROJECTED      } from '../../modules/local/agat_spstatistics'
include { AGAT_TO_MQC       as AGAT_PROJECTED_TO_MQC } from '../../modules/local/agat_to_mqc'

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

    // ── Statistics on projected models ───────────────────────────────────────
    // Pair each projected GFF with the (gunzipped) target genome FASTA so
    // AGAT can compute genome-coverage metrics via --gs. Join by target_id.
    ch_projected_agat_in = BAM_TO_GFF.out.gff3
        .map { meta, gff -> tuple(meta.target_id, meta + [kind: 'projected'], gff) }
        .combine(
            GUNZIP_TARGET.out.gunzip.map { meta, fa -> tuple(meta.id, fa) },
            by: 0
        )
        .map { _id, meta, gff, fa -> tuple(meta, gff, fa) }

    AGAT_PROJECTED(ch_projected_agat_in)
    AGAT_PROJECTED_TO_MQC(AGAT_PROJECTED.out.stats_yaml)

    // Keep SAMTOOLS_STATS.out.stats in the mix so MultiQC's native samtools
    // module still runs (it renders the "Percent mapped" bar chart). Our
    // curated 7-column TSV is rendered as a custom section nested under
    // Samtools (see assets/multiqc/sections.yml: parent_id: samtools).
    ch_mqc_files = AGAT_PROJECTED_TO_MQC.out.tsv
        .mix(SAMTOOLS_STATS.out.stats)
        .mix(SAMTOOLS_TO_MQC.out.tsv)

    emit:
    bam          = MINIMAP2_ALIGN.out.bam      // [ meta, *.bam     ]
    index        = MINIMAP2_ALIGN.out.index    // [ meta, *.bam.bai ]
    gff3         = BAM_TO_GFF.out.gff3         // [ meta, *.gff3    ]
    target_fasta = GUNZIP_TARGET.out.gunzip    // [ meta, fasta     ] × 1
    mqc_files    = ch_mqc_files                // [ meta, *_mqc.tsv ] × 16
}
