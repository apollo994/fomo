include { GUNZIP as GUNZIP_TARGET } from '../../modules/nf-core/gunzip/main'
include { MINIMAP2_INDEX          } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN          } from '../../modules/nf-core/minimap2/align/main'
include { BAM_TO_GFF              } from '../../modules/local/bam_to_gff'

workflow PROJECTION {
    take:
    ch_target            // [ meta(role:'target'), fasta.gz, gff3.gz ] × 1
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

    BAM_TO_GFF(MINIMAP2_ALIGN.out.bam)

    emit:
    bam   = MINIMAP2_ALIGN.out.bam    // [ meta, *.bam     ]
    index = MINIMAP2_ALIGN.out.index  // [ meta, *.bam.bai ]
    gff3  = BAM_TO_GFF.out.gff3       // [ meta, *.gff3    ]
}
