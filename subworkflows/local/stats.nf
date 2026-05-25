include { GUNZIP as GUNZIP_RAW_GFF } from '../../modules/nf-core/gunzip/main'
include { AGAT_SPSTATISTICS        } from '../../modules/local/agat_spstatistics'
include { AGAT_TO_MQC              } from '../../modules/local/agat_to_mqc'
include { SEQKIT_STATS             } from '../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_TO_MQC            } from '../../modules/local/seqkit_to_mqc'
include { MULTIQC                  } from '../../modules/nf-core/multiqc/main'

workflow STATS {
    take:
    ch_raw_source          // [ meta, fasta.gz, gff3.gz ] × N
    ch_raw_target          // [ meta, fasta.gz, gff3.gz ] × 1
    ch_filtered_gff3       // [ meta(decoy:false, feature_type), gff3 ] × 2N
    ch_decoy_gff3          // [ meta(decoy:true,  feature_type), gff3 ] × 2N
    ch_spliced_fasta       // [ meta(decoy:false, feature_type), fasta ] × 2N
    ch_decoy_spliced_fasta // [ meta(decoy:true,  feature_type), fasta ] × 2N
    ch_gffcompare_stats    // [ meta, *.stats ] × 32  (16 with-M + 16 no-M)

    main:

    // ── Decompress raw GFFs (source + target) ────────────────────────────────
    // AGAT does not handle .gff3.gz transparently — gunzip first.
    ch_raw_gff_gz = ch_raw_source.mix(ch_raw_target)
        .map { meta, fasta, gff3 -> tuple(meta + [kind: 'raw'], gff3) }

    GUNZIP_RAW_GFF(ch_raw_gff_gz)

    // ── AGAT statistics on raw + filtered + decoy GFFs ───────────────────────
    ch_agat_in = GUNZIP_RAW_GFF.out.gunzip
        .mix(ch_filtered_gff3.map { m, g -> tuple(m + [kind: 'filtered'], g) })
        .mix(ch_decoy_gff3   .map { m, g -> tuple(m + [kind: 'decoy'],    g) })

    AGAT_SPSTATISTICS(ch_agat_in)

    // ── Parse AGAT YAML into MultiQC custom-content TSVs ─────────────────────
    AGAT_TO_MQC(AGAT_SPSTATISTICS.out.stats_yaml)

    // ── SeqKit stats on every spliced FASTA (source + decoy) ─────────────────
    // FastQC was the user's first preference but FastQC 0.12.x does not
    // support FASTA inputs (defaults to FastQ for any file not ending in
    // .bam/.sam and fails on a `>`-prefixed header). SeqKit stats is the
    // standard nf-core tool for FASTA QC and MultiQC has native support.
    SEQKIT_STATS(ch_spliced_fasta.mix(ch_decoy_spliced_fasta))

    // Rewrite seqkit TSV first column so MultiQC's seqkit_stats custom_data
    // table uses a clean sample identifier (matching AGAT row names) instead
    // of the long FASTA filename.
    SEQKIT_TO_MQC(SEQKIT_STATS.out.stats)

    // ── MultiQC — assemble single tuple input per module signature ───────────
    // All three custom_data sections are ingested as TSVs by MultiQC config:
    //   • agat_summary  — curated 8-col headline table (one row per sample)
    //   • agat_full     — exhaustive per-feature-type metrics
    //   • seqkit_stats  — FASTA length/composition stats
    ch_mqc_files = SEQKIT_TO_MQC.out.mqc.map { _m, t -> t }
        .mix(AGAT_TO_MQC.out.summary.map { _m, t -> t })
        .mix(AGAT_TO_MQC.out.full   .map { _m, t -> t })
        .mix(ch_gffcompare_stats     .map { _m, s -> s })
        .collect()

    ch_mqc_input = ch_mqc_files.map { files ->
        tuple(
            [id: 'multiqc'],
            files,
            file(params.multiqc_config, checkIfExists: true),
            [],
            [],
            []
        )
    }

    MULTIQC(ch_mqc_input)

    emit:
    agat_stats_txt    = AGAT_SPSTATISTICS.out.stats_txt   // [ meta, *.stats.txt          ]
    agat_stats_yaml   = AGAT_SPSTATISTICS.out.stats_yaml  // [ meta, *.stats.yaml         ]
    agat_summary_tsv  = AGAT_TO_MQC.out.summary           // [ meta, *_summary_mqc.tsv    ]
    agat_full_tsv     = AGAT_TO_MQC.out.full              // [ meta, *_full_mqc.tsv       ]
    seqkit_stats      = SEQKIT_STATS.out.stats            // [ meta, *.spliced.tsv        ]
    multiqc_report    = MULTIQC.out.report                // [ meta, *.html               ]
    multiqc_data      = MULTIQC.out.data                  // [ meta, *_data               ]
}
