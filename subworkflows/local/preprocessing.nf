include { FILTER_TRANSCRIPT             } from '../../modules/local/filter_transcript'
include { GUNZIP as GUNZIP_FASTA        } from '../../modules/nf-core/gunzip/main'
include { GFFREAD as EXTRACT_SEQUENCES  } from '../../modules/nf-core/gffread/main'

workflow PREPROCESSING {
    take:
    ch_sources   // [ meta, fasta, gff3 ]

    main:

    // Decompress each source FASTA once (gffread needs uncompressed genome FASTA)
    GUNZIP_FASTA(
        ch_sources.map { meta, fasta, gff3 -> tuple(meta, fasta) }
    )

    // Filter each source GFF3 per feature_type (lnc_RNA, mRNA)
    ch_filter_input = ch_sources
        .combine(Channel.of('lnc_RNA', 'mRNA'))
        .map { meta, fasta, gff3, ftype ->
            tuple(meta + [feature_type: ftype], gff3, ftype)
        }

    FILTER_TRANSCRIPT(ch_filter_input)

    // Join filtered GFFs with their decompressed FASTA by species id (meta.id).
    // combine(by:0) acts as an inner join on the key element; multiMap keeps the
    // two outputs paired emission-by-emission for the GFFREAD signature.
    FILTER_TRANSCRIPT.out.gff3
        .map { meta, gff3 -> tuple(meta.id, meta, gff3) }
        .combine(
            GUNZIP_FASTA.out.gunzip.map { meta, fa -> tuple(meta.id, fa) },
            by: 0
        )
        .multiMap { id, meta, gff3, fa ->
            gff:   tuple(meta, gff3)
            fasta: fa
        }
        .set { ch_extract }

    EXTRACT_SEQUENCES(ch_extract.gff, ch_extract.fasta)

    emit:
    spliced_fasta = EXTRACT_SEQUENCES.out.gffread_fasta   // [ meta, fasta ]
    filtered_gff3 = FILTER_TRANSCRIPT.out.gff3            // [ meta, gff3  ]
}
