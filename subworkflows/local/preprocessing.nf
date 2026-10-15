include { FILTER_TRANSCRIPT                  } from '../../modules/local/filter_transcript'
include { GFF_TO_GENE_BED                    } from '../../modules/local/gff_to_gene_bed'
include { RELOCATE_LOCI                      } from '../../modules/local/relocate_loci'
include { RENAME_FASTA_HEADERS               } from '../../modules/local/rename_fasta_headers'
include { TD2_NONCODING                      } from './td2_noncoding'
include { FILTER_ANNOTATION                  } from '../../modules/local/filter_annotation'
include { GFFREAD as EXTRACT_SEQUENCES       } from '../../modules/nf-core/gffread/main'
include { GFFREAD as EXTRACT_DECOY_SEQUENCES } from '../../modules/nf-core/gffread/main'
include { FASTA_SIZES                        } from '../../modules/local/fasta_sizes'
include { BEDTOOLS_COMPLEMENT                } from '../../modules/nf-core/bedtools/complement/main'
include { GFF_STATS                          } from '../../modules/local/gff_stats'
include { GFF_STATS_TO_MQC                   } from '../../modules/local/gff_stats_to_mqc'
include { SEQKIT_STATS                       } from '../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_TO_MQC                      } from '../../modules/local/seqkit_to_mqc'

// Per source and target-independent: S tasks per step, never S·T.
workflow PREPROCESSING {
    take:
    ch_sources        // [ meta, fasta, gff3 ]
    feature_types     // List<String> of enabled feature types, derived once in workflows/fomo.nf.
                      // A plain List, expanded in a flatMap below (a Channel.value of a List
                      // would be spread by combine).

    main:

    // Inputs stay compressed (plans/26): each gffread task stages its genome on node-local
    // scratch, and every GFF3 consumer reads .gz. No samtools faidx: it rejects plain gzip.
    ch_source_fasta = ch_sources.map { meta, fasta, _gff3 -> tuple(meta.id, fasta) }

    // Candidates per (source, feature type): stranded, spliced and, with --longest_isoform,
    // the longest isoform per gene. They feed extraction, decoy relocation and
    // FILTER_ANNOTATION.
    ch_filter_input = ch_sources
        .flatMap { meta, _fasta, gff3 ->
            feature_types.collect { ft -> tuple(meta + [feature_type: ft, decoy: false], gff3, ft) }
        }

    FILTER_TRANSCRIPT(ch_filter_input)

    FILTER_TRANSCRIPT.out.gff3
        .map { meta, gff3 -> tuple(meta.id, meta, gff3) }
        .combine(ch_source_fasta, by: 0)
        .multiMap { _id, meta, gff3, fa ->
            gff:   tuple(meta, gff3)
            fasta: fa
        }
        .set { ch_extract }

    EXTRACT_SEQUENCES(ch_extract.gff, ch_extract.fasta)

    // Decoys (--include_decoy): source loci moved into the source's own intergenic space, the
    // false-positive baseline. Intergenic space is the complement of the RAW annotation's
    // gene-level features (coding genes included), whatever tracks are enabled, so the chain
    // runs once per source and its BED serves every feature type.
    if (params.include_decoy) {
        FASTA_SIZES(
            ch_sources.map { meta, fasta, _gff3 -> tuple(meta, fasta) }
        )

        // Gene BED in the chromosome order of the sizes file, as bedtools complement requires.
        ch_sources
            .map { meta, _fasta, gff3 -> tuple(meta.id, meta, gff3) }
            .combine(
                FASTA_SIZES.out.sizes.map { meta, sizes -> tuple(meta.id, sizes) },
                by: 0
            )
            .map { _id, meta, gff3, sizes -> tuple(meta, gff3, sizes) }
            .set { ch_gff_to_bed }

        GFF_TO_GENE_BED(ch_gff_to_bed)

        GFF_TO_GENE_BED.out.bed
            .map { meta, bed -> tuple(meta.id, meta, bed) }
            .combine(
                FASTA_SIZES.out.sizes.map { meta, sizes -> tuple(meta.id, sizes) },
                by: 0
            )
            .multiMap { _id, meta, bed, sizes ->
                bed:   tuple(meta, bed)
                sizes: sizes
            }
            .set { ch_complement }

        BEDTOOLS_COMPLEMENT(ch_complement.bed, ch_complement.sizes)

        FILTER_TRANSCRIPT.out.gff3
            .map { meta, gff3 -> tuple(meta.id, meta, gff3) }
            .combine(
                BEDTOOLS_COMPLEMENT.out.bed.map { meta, bed -> tuple(meta.id, bed) },
                by: 0
            )
            .map { _id, meta, gff3, bed -> tuple(meta + [decoy: true], gff3, bed) }
            .set { ch_relocate }

        RELOCATE_LOCI(ch_relocate)

        RELOCATE_LOCI.out.gff3
            .map { meta, gff3 -> tuple(meta.id, meta, gff3) }
            .combine(ch_source_fasta, by: 0)
            .multiMap { _id, meta, gff3, fa ->
                gff:   tuple(meta, gff3)
                fasta: fa
            }
            .set { ch_decoy_extract }

        EXTRACT_DECOY_SEQUENCES(ch_decoy_extract.gff, ch_decoy_extract.fasta)

        ch_decoy_gff3      = RELOCATE_LOCI.out.gff3
        ch_decoy_fasta_raw = EXTRACT_DECOY_SEQUENCES.out.gffread_fasta
    }
    else {
        ch_decoy_gff3      = channel.empty()
        ch_decoy_fasta_raw = channel.empty()
    }

    // `>tid|type|species` headers, for real and decoy tracks alike.
    RENAME_FASTA_HEADERS(
        EXTRACT_SEQUENCES.out.gffread_fasta
            .mix(ch_decoy_fasta_raw)
    )

    // TD2 drops real lncRNA with a complete, PSAURON-confident ORF; mRNA and decoys pass.
    RENAME_FASTA_HEADERS.out.fasta
        .branch { meta, _fa ->
            lnc:  !meta.decoy && meta.feature_type == 'lnc_RNA'
            rest: true
        }
        .set { ch_renamed }

    TD2_NONCODING(ch_renamed.lnc)

    ch_renamed_final = TD2_NONCODING.out.kept_fasta.mix(ch_renamed.rest)

    // The one filtered annotation per (source, feature type) is built FROM the FASTA that is
    // projected: FILTER_ANNOTATION keeps exactly the candidates whose renamed record survived,
    // matched through the id map by exact id, and fails if the three disagree. It is the
    // input "filtered" stats, the benchmarking reference of a 'both' species, and published.
    // Decoys come from the pre-TD2 candidates: they borrow only a model's shape.
    ch_final_fasta = ch_renamed_final.filter { meta, _fa -> !meta.decoy }

    // failOnMismatch: a source whose id map or FASTA is missing must abort the run.
    ch_annotation_in = FILTER_TRANSCRIPT.out.gff3
        .map { meta, gff -> tuple(meta.id, meta.feature_type, meta, gff) }
        .join(
            RENAME_FASTA_HEADERS.out.id_map
                .filter { meta, _map -> !meta.decoy }
                .map    { meta, map  -> tuple(meta.id, meta.feature_type, map) },
            by: [0, 1], failOnMismatch: true
        )
        .join(
            ch_final_fasta.map { meta, fa -> tuple(meta.id, meta.feature_type, fa) },
            by: [0, 1], failOnMismatch: true
        )
        .map { _id, _ft, meta, gff, map, fa -> tuple(meta, gff, map, fa) }

    FILTER_ANNOTATION(ch_annotation_in)

    ch_filtered_gff3 = FILTER_ANNOTATION.out.gff3

    // Stats of the raw, filtered and decoy annotations of each source (gff-feature-stats
    // reads .gz).
    ch_stats_gff = ch_sources.map { meta, _fa, gff3 -> tuple(meta + [kind: 'raw'], gff3) }
        .mix(ch_filtered_gff3.map { m, g -> tuple(m + [kind: 'filtered'], g) })
        .mix(ch_decoy_gff3          .map { m, g -> tuple(m + [kind: 'decoy'],    g) })

    GFF_STATS(ch_stats_gff)
    GFF_STATS_TO_MQC(GFF_STATS.out.json)

    SEQKIT_STATS(ch_renamed_final)
    SEQKIT_TO_MQC(SEQKIT_STATS.out.stats)

    ch_mqc_files = GFF_STATS_TO_MQC.out.tsv
        .mix(SEQKIT_TO_MQC.out.mqc)
        .mix(TD2_NONCODING.out.mqc)
        .mix(FILTER_TRANSCRIPT.out.mqc)

    // S sources, F feature types; the decoy emit is empty without --include_decoy.
    emit:
    spliced_fasta       = ch_final_fasta                                        // [ meta, fasta ] × F·S
    filtered_gff3       = ch_filtered_gff3                                      // [ meta, gff3  ] × F·S
    decoy_spliced_fasta = ch_renamed_final.filter { meta, _fa -> meta.decoy }   // [ meta, fasta ] × F·S
    mqc_files           = ch_mqc_files                                          // [ meta, path  ]
}
