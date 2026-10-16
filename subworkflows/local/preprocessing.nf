include { FILTER_TRANSCRIPT                  } from '../../modules/local/filter_transcript'
include { GFF_TO_GENE_BED                    } from '../../modules/local/gff_to_gene_bed'
include { RELOCATE_LOCI                      } from '../../modules/local/relocate_loci'
include { RENAME_FASTA_HEADERS               } from '../../modules/local/rename_fasta_headers'
include { TD2_NONCODING                      } from './td2_noncoding'
include { FILTER_GFF_BY_ID as FILTER_ANNOTATION         } from '../../modules/local/filter_gff_by_id'
include { EXTRACT_TRANSCRIPTS as EXTRACT_SEQUENCES       } from '../../modules/local/extract_transcripts'
include { EXTRACT_TRANSCRIPTS as EXTRACT_DECOY_SEQUENCES } from '../../modules/local/extract_transcripts'
include { FASTA_SIZES                        } from '../../modules/local/fasta_sizes'
include { BEDTOOLS_COMPLEMENT                } from '../../modules/local/bedtools_complement'
include { GFF_STATS                          } from '../../modules/local/gff_stats'
include { SEQKIT_STATS                       } from '../../modules/local/seqkit_stats'
include { withTrack; pairById                } from './fomo_utils'

// Per source and target-independent: S tasks per step, never S·T.
workflow PREPROCESSING {
    take:
    ch_sources        // [ meta, fasta, gff3 ]
    feature_types     // List<String> of enabled feature types, derived once in workflows/fomo.nf.
                      // A plain List, expanded in a flatMap below (a Channel.value of a List
                      // would be spread by combine).

    main:

    // Inputs stay compressed: EXTRACT_TRANSCRIPTS stages its genome on node-local scratch, and
    // every GFF3 consumer reads .gz. No samtools faidx: it rejects plain gzip. Rationale: plans/26.
    ch_source_fasta = ch_sources.map { meta, fasta, _gff3 -> tuple(meta.id, fasta) }

    // Candidates per (source, track): stranded, spliced and, with --longest_isoform, the longest
    // isoform per gene. They feed extraction, decoy relocation and FILTER_ANNOTATION.
    ch_filter_input = ch_sources
        .flatMap { meta, _fasta, gff3 ->
            feature_types.collect { ft -> tuple(withTrack(meta + [feature_type: ft, decoy: false]), gff3) }
        }

    FILTER_TRANSCRIPT(ch_filter_input)

    EXTRACT_SEQUENCES(pairById(FILTER_TRANSCRIPT.out.gff3, ch_source_fasta))

    // Decoys (--include_decoy): source loci moved into the source's own intergenic space, the
    // false-positive baseline. Intergenic space is the complement of the RAW annotation's
    // gene-level features (coding genes included), whatever tracks are enabled, so the chain
    // runs once per source and its BED serves every feature type. Without the flag the chain
    // gets no input and runs no task.
    ch_decoy_sources = params.include_decoy ? ch_sources : channel.empty()

    FASTA_SIZES(ch_decoy_sources.map { meta, fasta, _gff3 -> tuple(meta, fasta) })

    ch_sizes = FASTA_SIZES.out.sizes.map { meta, sizes -> tuple(meta.id, sizes) }

    // Gene BED in the chromosome order of the sizes file, as bedtools complement requires.
    GFF_TO_GENE_BED(pairById(ch_decoy_sources.map { meta, _fasta, gff3 -> tuple(meta, gff3) }, ch_sizes))

    BEDTOOLS_COMPLEMENT(pairById(GFF_TO_GENE_BED.out.bed, ch_sizes))

    RELOCATE_LOCI(
        pairById(
            FILTER_TRANSCRIPT.out.gff3.map { meta, gff3 -> tuple(withTrack(meta + [decoy: true]), gff3) },
            BEDTOOLS_COMPLEMENT.out.bed.map { meta, bed -> tuple(meta.id, bed) }
        )
    )

    EXTRACT_DECOY_SEQUENCES(pairById(RELOCATE_LOCI.out.gff3, ch_source_fasta))

    // `>tid|type|species` headers, for real and decoy tracks alike.
    RENAME_FASTA_HEADERS(EXTRACT_SEQUENCES.out.fasta.mix(EXTRACT_DECOY_SEQUENCES.out.fasta))

    // TD2 drops real lncRNA with a complete, PSAURON-confident ORF; mRNA and decoys pass.
    RENAME_FASTA_HEADERS.out.fasta
        .branch { meta, _fa ->
            lnc:  meta.td2
            rest: true
        }
        .set { ch_renamed }

    TD2_NONCODING(ch_renamed.lnc.map { meta, fa -> tuple(meta + [stage: 'input'], fa) })

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
        .mix(ch_filtered_gff3.map       { m, g -> tuple(m + [kind: 'filtered'], g) })
        .mix(RELOCATE_LOCI.out.gff3.map { m, g -> tuple(m + [kind: 'decoy'],    g) })

    GFF_STATS(ch_stats_gff.map { meta, gff3 -> tuple(meta, gff3, []) })

    SEQKIT_STATS(ch_renamed_final)

    ch_data = FILTER_TRANSCRIPT.out.data
        .mix(GFF_STATS.out.data)
        .mix(SEQKIT_STATS.out.data)
        .mix(TD2_NONCODING.out.data)

    // S sources, F feature types; the decoy emit is empty without --include_decoy.
    emit:
    spliced_fasta       = ch_final_fasta                                        // [ meta, fasta ] × F·S
    filtered_gff3       = ch_filtered_gff3                                      // [ meta, gff3  ] × F·S
    decoy_spliced_fasta = ch_renamed_final.filter { meta, _fa -> meta.decoy }   // [ meta, fasta ] × F·S
    data                = ch_data                                               // [ meta, *.data.json ]
}
