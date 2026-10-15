process GFF_STATS {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_single'
    label 'img_ubuntu'
    label 'amd64_only'

    // bin/gff-feature-stats (provenance: bin/gff-feature-stats.PROVENANCE) reads .gff3.gz; its
    // input must be grouped by seqid.
    input:
    tuple val(meta), path(gff)          // .gff3 or .gff3.gz; no genome FASTA needed

    output:
    tuple val(meta), path("*.stats.json"), emit: json
    tuple val("${task.process}"), val('gff-feature-stats'), eval("gff-feature-stats -V | sed 's/gff-feature-stats //'"), topic: versions, emit: versions_gff_feature_stats

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
