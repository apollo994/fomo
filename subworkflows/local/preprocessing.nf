include { FILTER_TRANSCRIPT                  } from '../../modules/local/filter_transcript'
include { GFF_TO_GENE_BED                    } from '../../modules/local/gff_to_gene_bed'
include { RELOCATE_LOCI                      } from '../../modules/local/relocate_loci'
include { RENAME_FASTA_HEADERS               } from '../../modules/local/rename_fasta_headers'
include { GUNZIP as GUNZIP_FASTA             } from '../../modules/nf-core/gunzip/main'
include { GFFREAD as EXTRACT_SEQUENCES       } from '../../modules/nf-core/gffread/main'
include { GFFREAD as EXTRACT_DECOY_SEQUENCES } from '../../modules/nf-core/gffread/main'
include { SAMTOOLS_FAIDX                     } from '../../modules/nf-core/samtools/faidx/main'
include { BEDTOOLS_COMPLEMENT                } from '../../modules/nf-core/bedtools/complement/main'

workflow PREPROCESSING {
    take:
    ch_sources   // [ meta, fasta, gff3 ]

    main:

    // ── Source path ──────────────────────────────────────────────────────────
    // Decompress each source FASTA once (gffread needs uncompressed genome)
    GUNZIP_FASTA(
        ch_sources.map { meta, fasta, gff3 -> tuple(meta, fasta) }
    )

    // Filter each source GFF3 per feature_type (lnc_RNA, mRNA)
    ch_filter_input = ch_sources
        .combine(Channel.of('lnc_RNA', 'mRNA'))
        .map { meta, fasta, gff3, ftype ->
            tuple(meta + [feature_type: ftype, decoy: false], gff3, ftype)
        }

    FILTER_TRANSCRIPT(ch_filter_input)

    // Join filtered GFFs with their decompressed FASTA by species id
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

    // ── Per-source intergenic intervals ──────────────────────────────────────
    // Index each decompressed source FASTA to get chromosome sizes
    SAMTOOLS_FAIDX(
        GUNZIP_FASTA.out.gunzip.map { meta, fa -> tuple(meta, fa, []) },
        true
    )

    // Extract gene-level features from each source GFF3 → sorted BED
    GFF_TO_GENE_BED(
        ch_sources.map { meta, fasta, gff3 -> tuple(meta, gff3) }
    )

    // Pair each source gene BED with its chromosome sizes, then complement
    GFF_TO_GENE_BED.out.bed
        .map { meta, bed -> tuple(meta.id, meta, bed) }
        .combine(
            SAMTOOLS_FAIDX.out.sizes.map { meta, sizes -> tuple(meta.id, sizes) },
            by: 0
        )
        .multiMap { id, meta, bed, sizes ->
            bed:   tuple(meta, bed)
            sizes: sizes
        }
        .set { ch_complement }

    BEDTOOLS_COMPLEMENT(ch_complement.bed, ch_complement.sizes)

    // ── Decoy path ───────────────────────────────────────────────────────────
    // Pair each filtered source GFF3 with its source-species intergenic BED
    FILTER_TRANSCRIPT.out.gff3
        .map { meta, gff3 -> tuple(meta.id, meta, gff3) }
        .combine(
            BEDTOOLS_COMPLEMENT.out.bed.map { meta, bed -> tuple(meta.id, bed) },
            by: 0
        )
        .map { id, meta, gff3, bed -> tuple(meta + [decoy: true], gff3, bed) }
        .set { ch_relocate }

    RELOCATE_LOCI(ch_relocate)

    // Extract spliced decoy sequences from each source FASTA
    RELOCATE_LOCI.out.gff3
        .map { meta, gff3 -> tuple(meta.id, meta, gff3) }
        .combine(
            GUNZIP_FASTA.out.gunzip.map { meta, fa -> tuple(meta.id, fa) },
            by: 0
        )
        .multiMap { id, meta, gff3, fa ->
            gff:   tuple(meta, gff3)
            fasta: fa
        }
        .set { ch_decoy_extract }

    EXTRACT_DECOY_SEQUENCES(ch_decoy_extract.gff, ch_decoy_extract.fasta)

    // ── Rename FASTA headers ─────────────────────────────────────────────────
    // Mix source and decoy FASTAs; meta carries feature_type and decoy flag so
    // the module derives the correct type label without separate invocations.
    RENAME_FASTA_HEADERS(
        EXTRACT_SEQUENCES.out.gffread_fasta
            .mix(EXTRACT_DECOY_SEQUENCES.out.gffread_fasta)
    )

    emit:
    spliced_fasta       = RENAME_FASTA_HEADERS.out.fasta.filter { meta, fa -> !meta.decoy } // [ meta, fasta ] × 2N
    filtered_gff3       = FILTER_TRANSCRIPT.out.gff3                                        // [ meta, gff3  ] × 2N
    intergenic_bed      = BEDTOOLS_COMPLEMENT.out.bed                                       // [ meta, bed   ] × N
    decoy_gff3          = RELOCATE_LOCI.out.gff3                                            // [ meta, gff3  ] × 2N
    decoy_spliced_fasta = RENAME_FASTA_HEADERS.out.fasta.filter { meta, fa -> meta.decoy }  // [ meta, fasta ] × 2N
}
