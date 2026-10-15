process GFF_STATS {
    tag "${meta.id}.${meta.kind}"
    label 'process_single'

    // bin/gff-feature-stats (vendored from github.com/apollo994/gff-feature-stats @ 67ffea7,
    // reports 0.2.0) reads .gff3.gz. Inputs
    // must be grouped by seqid (enforced above 200k lines). The binary is x86-64 and needs
    // glibc >= 2.34 (ubuntu:22.04 has 2.35); the docker profile pins amd64 for it.
    input:
    tuple val(meta), path(gff)          // .gff3 or .gff3.gz; no genome FASTA needed

    output:
    tuple val(meta), path("*.stats.json"), emit: json
    tuple val("${task.process}"), val('gff-feature-stats'), eval("gff-feature-stats -V | sed 's/gff-feature-stats //'"), topic: versions, emit: versions_gff_feature_stats

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    // A feature-less GFF yields empty stats and exit 0.
    """
    gff-feature-stats ${gff} ${prefix}.stats.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    """
    echo '{"gene_category_stats":{},"transcript_type_stats":{}}' > ${prefix}.stats.json
    """
}
