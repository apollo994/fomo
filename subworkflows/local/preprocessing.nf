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

workflow PREPROCESSING {
    take:
    ch_sources        // [ meta, fasta, gff3 ]
    feature_types     // plain List<String> of enabled feature types (lnc_RNA [, mRNA]).
                      // Derived ONCE in workflows/fomo.nf and shared with BENCHMARKING —
                      // never re-derive it here (see the comment at that call site).
                      // A List, not a channel, on purpose: it is expanded in a flatMap
                      // closure below. See the fomo.nf comment for why .combine() on a
                      // Channel.value(List) is wrong here.

    main:

    // ── Source path ──────────────────────────────────────────────────────────
    // The source FASTA stays as the samplesheet gives it (normally plain gzip) — there
    // is no GUNZIP task (plans/26). gffread's `-g` reader has no gzip path, so each
    // gffread task (EXTRACT_SEQUENCES / EXTRACT_DECOY_SEQUENCES) decompresses its genome
    // into node-local scratch itself (FOMO-local branch of modules/nf-core/gffread):
    // a shared plain copy on NFS was read twice per task and stalled large assemblies.
    // The decoy chain's chromosome sizes come from FASTA_SIZES, one streaming pass.
    // Do NOT reach for samtools faidx on the .gz: it accepts bgzip but hard-errors on
    // plain gzip, and while the committed test FASTAs are bgzip
    // (bin/subsample_test_data.sh), the real Ensembl/NCBI assemblies are plain gzip —
    // such a shortcut passes -profile test and fails on the cluster.
    // The GFF3 side reads .gz too: FILTER_TRANSCRIPT and GFF_TO_GENE_BED decompress
    // inline (to local scratch / a pipe), gff-feature-stats natively.
    ch_source_fasta = ch_sources.map { meta, fasta, _gff3 -> tuple(meta.id, fasta) }

    // Filter each source GFF3 per enabled feature_type (lnc_RNA [, mRNA]): spliced
    // transcripts only and, with --longest_isoform (ext.args), one per gene — the
    // longest. Its output is the CANDIDATE set: it feeds extraction (→ TD2 →
    // projection), decoy relocation, and FILTER_ANNOTATION, which turns it into the
    // one filtered annotation everything downstream reads.
    ch_filter_input = ch_sources
        .flatMap { meta, _fasta, gff3 ->
            feature_types.collect { ft -> tuple(meta + [feature_type: ft, decoy: false], gff3, ft) }
        }

    FILTER_TRANSCRIPT(ch_filter_input)

    // Join filtered GFFs with their (still compressed) source FASTA by species id
    FILTER_TRANSCRIPT.out.gff3
        .map { meta, gff3 -> tuple(meta.id, meta, gff3) }
        .combine(ch_source_fasta, by: 0)
        .multiMap { _id, meta, gff3, fa ->
            gff:   tuple(meta, gff3)
            fasta: fa
        }
        .set { ch_extract }

    EXTRACT_SEQUENCES(ch_extract.gff, ch_extract.fasta)

    // ── Decoy track (opt-in: --include_decoy) ────────────────────────────────
    // Decoys are the false-positive baseline: source loci relocated into that
    // source's own intergenic space, so anything that still projects onto the
    // target is a false positive. The intergenic machinery (sizes → gene BED →
    // complement) exists only to feed RELOCATE_LOCI, so the whole branch is
    // gated together. One decoy track is built per enabled feature_type.
    //
    // INVARIANT — decoy DESTINATIONS ignore the feature-type list. GFF_TO_GENE_BED
    // is fed from ch_sources (the RAW samplesheet GFF3) and keeps every gene-level
    // feature (gene | pseudogene | ncRNA_gene), so intergenic space is the
    // complement of the ENTIRE annotation, coding genes included, whether or not
    // mRNA is being transferred. The chain is per SOURCE (N tasks, not F·N) and the
    // one BED per species is broadcast to every feature_type below. So do NOT gate
    // any of it on params.include_mrna, and do NOT "simplify" GFF_TO_GENE_BED's
    // input to the per-feature-type FILTER_TRANSCRIPT output — that would redefine
    // intergenic as "outside lncRNA genes" and let lncRNA decoys land inside real
    // mRNA loci.
    if (params.include_decoy) {
        // Chromosome sizes of each source genome, streamed from the .gz (plans/26)
        FASTA_SIZES(
            ch_sources.map { meta, fasta, _gff3 -> tuple(meta, fasta) }
        )

        // Extract gene-level features from each source GFF3 → BED, sorted in the
        // chromosome order of the FASTA-derived sizes file (required so bedtools
        // complement downstream accepts the input).
        ch_sources
            .map { meta, _fasta, gff3 -> tuple(meta.id, meta, gff3) }
            .combine(
                FASTA_SIZES.out.sizes.map { meta, sizes -> tuple(meta.id, sizes) },
                by: 0
            )
            .map { _id, meta, gff3, sizes -> tuple(meta, gff3, sizes) }
            .set { ch_gff_to_bed }

        GFF_TO_GENE_BED(ch_gff_to_bed)

        // Pair each source gene BED with its chromosome sizes, then complement
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

        // Pair each filtered source GFF3 with its source-species intergenic BED
        FILTER_TRANSCRIPT.out.gff3
            .map { meta, gff3 -> tuple(meta.id, meta, gff3) }
            .combine(
                BEDTOOLS_COMPLEMENT.out.bed.map { meta, bed -> tuple(meta.id, bed) },
                by: 0
            )
            .map { _id, meta, gff3, bed -> tuple(meta + [decoy: true], gff3, bed) }
            .set { ch_relocate }

        RELOCATE_LOCI(ch_relocate)

        // Extract spliced decoy sequences from each source FASTA
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
        ch_intergenic_bed  = BEDTOOLS_COMPLEMENT.out.bed
    }
    else {
        ch_decoy_gff3      = channel.empty()
        ch_decoy_fasta_raw = channel.empty()
        ch_intergenic_bed  = channel.empty()
    }

    // ── Rename FASTA headers ─────────────────────────────────────────────────
    // Mix source and decoy FASTAs; meta carries feature_type and decoy flag so
    // the module derives the correct type label without separate invocations.
    RENAME_FASTA_HEADERS(
        EXTRACT_SEQUENCES.out.gffread_fasta
            .mix(ch_decoy_fasta_raw)
    )

    // ── TD2 coding-potential filter (real lncRNA only) ────────────────────────
    // Split the renamed FASTAs: real lncRNA go through TD2; everything else
    // (real mRNA + all decoys) passes through untouched. A real lncRNA with a
    // complete, PSAURON-confident ORF is dropped (coding potential).
    RENAME_FASTA_HEADERS.out.fasta
        .branch { meta, _fa ->
            lnc:  !meta.decoy && meta.feature_type == 'lnc_RNA'
            rest: true
        }
        .set { ch_renamed }

    TD2_NONCODING(ch_renamed.lnc)

    // Non-coding lncRNA ⊕ untouched (real mRNA + decoys) → the full renamed set
    // with lncRNA coding-filtered. spliced_fasta / decoy_spliced_fasta and
    // SeqKit stats are all derived from this so they reflect the filtering.
    ch_renamed_final = TD2_NONCODING.out.kept_fasta.mix(ch_renamed.rest)

    // ── The one filtered annotation per (source, feature_type) ──────────────
    // <id>.<ft>.filtered.gff3 (plans/20_longest_isoform.md) is built FROM the FASTA
    // that is projected, not alongside it: FILTER_ANNOTATION keeps exactly the
    // candidate transcripts whose renamed record survived into ch_final_fasta (TD2
    // non-coding for lncRNA, every record for mRNA), translated back through
    // RENAME_FASTA_HEADERS' id map by EXACT id. It fails loudly if the three disagree,
    // so this file and spliced_fasta describe the same set by construction. It is
    // then used three ways: input "filtered" stats here, the gffcompare reference in
    // BENCHMARKING (for 'both' species), and published under sources/<id>/annotation/.
    // Projection does not wait for it — spliced_fasta is emitted straight from
    // ch_final_fasta.
    //
    // Decoys deliberately do NOT come from here: RELOCATE_LOCI (above) reads the
    // FILTER_TRANSCRIPT candidates. A decoy borrows only a model's SHAPE (exon count,
    // exon and intron lengths) and takes its sequence from intergenic DNA at the new
    // site, so TD2's verdict on the original locus says nothing about it — and not
    // waiting on TD2 keeps the decoy branch off its critical path. Decoys do inherit
    // the longest-isoform selection, since that happens inside FILTER_TRANSCRIPT.
    ch_final_fasta = ch_renamed_final.filter { meta, _fa -> !meta.decoy }

    // Keyed on [id, feature_type]: one candidate GFF3, one id map and one final FASTA
    // per source and feature type. failOnMismatch — a source whose map or FASTA never
    // arrives must abort the run, not silently vanish from the results, the
    // benchmarking reference and the report.
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

    // ── Statistics & MultiQC adapters ─────────────────────────────────────────
    // Union of raw/filtered/decoy GFFs, all on the source genome. The raw source
    // GFF3 goes in still-gzipped — gff-feature-stats reads .gz directly. The
    // "filtered" set is the one filtered annotation (spliced, longest isoform, lncRNA
    // TD2-noncoding), so the input stats describe exactly what was projected.
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
        .mix(FILTER_TRANSCRIPT.out.mqc)   // transcript funnel: spliced / longest isoform

    // Cardinalities below use N = sources, F = enabled feature types (1, or 2 with
    // --include_mrna). The decoy emits are all EMPTY without --include_decoy.
    emit:
    spliced_fasta       = ch_final_fasta                                       // [ meta, fasta ] × F·N (lncRNA coding-filtered)
    filtered_gff3       = ch_filtered_gff3                                     // [ meta, gff3  ] × F·N — describes spliced_fasta exactly; the benchmarking reference
    intergenic_bed      = ch_intergenic_bed                                    // [ meta, bed   ] × N   (0 without --include_decoy)
    decoy_gff3          = ch_decoy_gff3                                        // [ meta, gff3  ] × F·N (0 without --include_decoy)
    decoy_spliced_fasta = ch_renamed_final.filter { meta, _fa -> meta.decoy }   // [ meta, fasta ] × F·N (0 without --include_decoy)
    mqc_files           = ch_mqc_files                                         // [ meta, path  ]
}
