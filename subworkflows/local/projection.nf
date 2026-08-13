include { MAYBE_GUNZIP as GUNZIP_TARGET            } from '../../modules/local/maybe_gunzip'
include { MERGE_FASTA as MERGE_SOURCE_FASTA        } from '../../modules/local/merge_fasta'
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
include { FILTER_GFF_BY_ID  as FILTER_ALLMODELS    } from '../../modules/local/filter_gff_by_id'
include { GFF_BY_SOURCE     as SPLIT_GFF_BY_SOURCE } from '../../modules/local/gff_by_source'

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

    // ── Merge every source's spliced transcripts, per feature class ───────────
    // One FASTA per gene type — {lnc_RNA, mRNA, lnc_RNA_decoy, mRNA_decoy}, of which
    // only lnc_RNA exists by default — so a target is aligned ONCE per class instead
    // of once per (source, class). No information is lost: every record already
    // carries its species in the header (`>tid|type|species`, RENAME_FASTA_HEADERS),
    // minimap2 copies that into the BAM QNAME, and bam_to_gff.sh turns it into a
    // `source=` GFF attribute — which is how the per-source view is recovered below.
    //
    // Classes stay separate rather than merging into one global FASTA because every
    // downstream grouping (consensus, benchmarking, MultiQC ordering) is keyed on
    // gene type; splitting them back out would need a second demux dimension.
    //
    // This is target-independent — F·D tasks, NOT F·D·T. Do not fan it out per target.
    ch_merged = ch_spliced_fasta
        .mix(ch_decoy_spliced)
        .map { meta, fa ->
            def gtype = meta.feature_type + (meta.decoy ? '_decoy' : '')
            tuple([id: gtype, feature_type: meta.feature_type, decoy: meta.decoy, gtype: gtype], fa)
        }
        .groupTuple()
        .map { meta, fas -> tuple(meta, fas.sort { it.name }) }   // deterministic concat order

    MERGE_SOURCE_FASTA(ch_merged)

    // Fan the merged tracks out over targets: F·D·T alignments. combine() without
    // `by` is a full cartesian product and buffers the right-hand side, so each of
    // the T references is reused across all F·D tracks. Self-projections are part of
    // the merged FASTA like any other source — they are the Sn/Pr ceiling control,
    // and they stay in allModels (they are only excluded from the top-N ranking pool,
    // in consensus_top.nf).
    //
    // multiMap, not a value channel: with T > 1 the reference must be paired with
    // its own reads task-by-task. A broadcast/.first() reference would silently
    // align every source onto whichever target's index materialised first.
    ch_align = MERGE_SOURCE_FASTA.out.fasta
        .combine(ch_ref)
        .multiMap { meta, reads, tid, mmi, _fa ->
            reads:     tuple(meta + [id: 'allModels', target_id: tid], reads)
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
    // plot and can't be customised via config). One row per (target, gene type)
    // — the merged BAM has no per-source breakdown; per-source model counts come
    // from GFF_STATS_PROJECTED instead.
    SAMTOOLS_TO_MQC(SAMTOOLS_STATS.out.stats)

    BAM_TO_GFF(MINIMAP2_ALIGN.out.bam)

    // ── TD2 coding-potential filter on the projected lncRNA ──────────────────
    // Runs ONCE per target (T tasks, not S·T) on the whole all-sources model set —
    // the split happens after it, so every per-source GFF3 and the consensus inherit
    // the same filtering. Real mRNA and all decoys pass through untouched (they still
    // go through FILTER_ALLMODELS with an empty drop list, so that a single process
    // produces every allModels.*.raw.gff3).
    BAM_TO_GFF.out.gff3
        .branch { meta, _g ->
            lnc:  !meta.decoy && meta.feature_type == 'lnc_RNA'
            rest: true
        }
        .set { ch_proj }

    // The transcript FASTA must be extracted from ITS OWN target genome — hence the
    // join on target_id rather than a broadcast FASTA. Getting this wrong is silent:
    // gffread would happily read the coordinates of target B out of target A's
    // sequence and hand TD2 nonsense to score.
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

    // lncRNA gets TD2's coding ids; every other class gets an empty drop list so it
    // travels the same path and comes out of the same process under the same naming
    // convention. filter_gff_by_id.py treats an empty file as "drop nothing".
    ch_filter_in = ch_proj.lnc
        .join(TD2_NONCODING_PROJ.out.coding_ids)
        .mix(
            ch_proj.rest.map { meta, gff ->
                tuple(meta, gff, file("${projectDir}/assets/empty_ids.txt", checkIfExists: true))
            }
        )

    FILTER_ALLMODELS(ch_filter_in)

    ch_allmodels_raw = FILTER_ALLMODELS.out.gff3   // [ meta(id:'allModels'), gff3 ] × F·D·T

    // ── Recover the per-source view ──────────────────────────────────────────
    // A pure attribute filter on `source=`, so each output is byte-for-byte the
    // subset of allModels belonging to that species — and, because the input is
    // coordinate-sorted, still grouped by seqid as gff-feature-stats requires.
    // [] = split mode (the third slot carries a keep-list in subset mode; see
    // consensus_top.nf).
    SPLIT_GFF_BY_SOURCE(
        ch_allmodels_raw.map { meta, gff -> tuple(meta, gff, []) }
    )

    // Re-key each split file back to its source species. The filename was built by
    // the same closure that supplies --name-template (conf/modules.config), so the
    // prefix and suffix are known exactly — strip them by length rather than pattern
    // matching, which would be fragile against species names full of '_' and digits.
    ch_projected_persource = SPLIT_GFF_BY_SOURCE.out.gff3
        .transpose()
        .map { meta, gff ->
            def pre = "${meta.target_id}.from_"
            def suf = ".${meta.feature_type}${meta.decoy ? '.decoy' : ''}.projected.gff3"
            tuple([id:           gff.name[pre.size()..-(suf.size() + 1)],
                   feature_type: meta.feature_type,
                   decoy:        meta.decoy,
                   target_id:    meta.target_id], gff)
        }

    // ── Collapse the raw model set into a consensus ───────────────────────────
    // gffcompare in combine mode (no reference annotation) over the SINGLE
    // allModels GFF3 — it merges transcripts sharing an intron chain, so one input
    // file holding every source collapses exactly as N per-source files used to.
    // One consensus per (target, gene type): <target>.allModels.<gtype>.combined.gtf.
    //
    // Self-projections are INCLUDED, unlike the pre-single-mapping consensus. The
    // rule now lives in exactly one place: self is part of allModels, and is
    // excluded only from the top-N RANKING pool (consensus_top.nf), where it would
    // otherwise take the #1 F1 slot on every annotated target and turn top3 into a
    // restatement of the target's own annotation. Do not re-add a self filter here.
    //
    // `gtype` rides in on the merge meta (see ch_merged above), so no re-derivation.
    GFFCOMPARE_COMBINE(
        ch_allmodels_raw,
        [[:], [], []],   // no reference sequence  (-s)
        [[:], []]        // no reference annotation (-r) → pure combine mode
    )

    // Convert each combined GTF → GFF3 so it matches the rest of the pipeline, and
    // give it the pseudo-source id 'allModels_collapsed' so it flows through the same
    // GFF stats / benchmarking naming closures as per-source projections (yielding
    // "from_allModels_collapsed" rows).
    COMBINED_GTF_TO_GFF(
        GFFCOMPARE_COMBINE.out.combined_gtf.map { meta, gtf -> tuple(meta + [id: 'allModels_collapsed'], gtf) },
        []   // no genome FASTA needed for a pure GTF→GFF3 conversion
    )

    // Union of per-source projections + the two aggregate models. Everything
    // downstream (GFF stats, benchmarking, the accuracy scatter) treats each
    // aggregate as one more pseudo-source, so raw and collapsed are scored
    // side by side against the per-source models.
    ch_projected = ch_projected_persource
        .mix(ch_allmodels_raw.map { meta, gff -> tuple(meta + [id: 'allModels_raw'], gff) })
        .mix(COMBINED_GTF_TO_GFF.out.gffread_gff)

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
    bam                 = MINIMAP2_ALIGN.out.bam     // [ meta, *.bam     ] × F·D·T (all-sources)
    index               = MINIMAP2_ALIGN.out.index   // [ meta, *.bam.bai ] × F·D·T
    allmodels_raw       = ch_allmodels_raw           // [ meta(id:'allModels'), *.gff3 ] × F·D·T (self included)
    allmodels_collapsed = COMBINED_GTF_TO_GFF.out.gffread_gff  // [ meta(id:'allModels_collapsed'), *.gff3 ] × F·D·T
    gff3                = ch_projected               // [ meta, *.gff3 ] × F·S·D·T + 2·F·D·T (per-source + both aggregates)
    mqc_files           = ch_mqc_files               // [ meta, *_mqc.tsv ] (one table per projected model + samtools + TD2)
}
