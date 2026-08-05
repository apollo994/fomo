process GFF_STATS {
    tag "${meta.id}.${meta.kind}"
    label 'process_single'

    // gff-feature-stats replaces agat_sp_statistics.pl: same gene-category and
    // transcript-type metrics plus intron statistics, at a flat ~20 MB RSS and
    // ~1 s for a 660 MB GFF (AGAT needed up to 12 GB / 2 h here). It also reads
    // .gff3.gz directly, so no gunzip step is needed upstream.
    //
    // The tool is vendored as bin/gff-feature-stats — a release build of
    // gff-feature-stats v0.2.0 (github.com/apollo994/gff-feature-stats @ 67ffea7).
    // Nextflow binds bin/ into the task container and prepends it to PATH.
    // NOTE: the commit is the real provenance — 67ffea7 did not bump the crate
    // version, so `-V` (and the versions topic) reports 0.2.0 for it and for the
    // earlier 6b31e3a build alike.
    //
    // 67ffea7 added an input requirement: each seqid's records must be
    // CONTIGUOUS (order within a seqid is free). Interleaved seqids are rejected
    // with exit 1, but only once the input crosses the tool's 200k-line batch
    // threshold — smaller files are a single batch and always accepted. Every
    // GFF fed to this process is grouped: Ensembl input, gffread output, and
    // BAM_TO_GFF output (MINIMAP2_ALIGN sorts via samtools). Decoys are grouped
    // explicitly by bin/relocate_loci.py, which otherwise interleaves them.
    //
    // The binary is dynamically linked and needs glibc >= 2.34 (ubuntu:22.04 ships
    // 2.35 — verified). Do NOT point `container` at an older base image without
    // re-testing, or rebuild statically:
    //   RUSTFLAGS='-C target-feature=+crt-static' cargo build --release
    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(gff)          // .gff3 or .gff3.gz; no genome FASTA needed

    output:
    tuple val(meta), path("*.stats.json"), emit: json
    tuple val("${task.process}"), val('gff-feature-stats'), eval("gff-feature-stats -V | sed 's/gff-feature-stats //'"), topic: versions, emit: versions_gff_feature_stats

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    // A feature-less GFF yields empty stats maps and exit 0, so no guard is
    // needed here (AGAT required one). Progress and any span-divergence
    // warnings go to stderr → .command.err.
    """
    gff-feature-stats ${gff} ${prefix}.stats.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    """
    echo '{"gene_category_stats":{},"transcript_type_stats":{}}' > ${prefix}.stats.json
    """
}
