process FILTER_TRANSCRIPT {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_single'
    label 'img_python'

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*.gff3")                        , emit: gff3
    // One-row funnel (n_transcripts, n_stranded, n_spliced, n_selected), Sample <id>.<ft>.
    tuple val(meta), path("*_transcript_filter_mqc.tsv")   , emit: mqc
    tuple val(meta), path("*.transcript_filter.data.json") , emit: data
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    def args   = task.ext.args ?: ''            // --longest (params.longest_isoform)
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    """
    filter_transcript.py \\
        ${args} \\
        --gff ${gff3} \\
        --feature-type ${meta.feature_type} \\
        --counts ${prefix}_transcript_filter_mqc.tsv \\
        --sample ${meta.id}.${meta.feature_type} \\
        --record ${prefix}.transcript_filter.data.json \\
        --identity '${task.ext.identity}' \\
        --output ${prefix}.gff3
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    """
    touch ${prefix}.gff3 ${prefix}.transcript_filter.data.json
    printf 'Sample\\tfeature\\tn_transcripts\\tn_stranded\\tn_spliced\\tn_selected\\n${meta.id}.${meta.feature_type}\\t${meta.feature_type}\\t0\\t0\\t0\\t0\\n' > ${prefix}_transcript_filter_mqc.tsv
    """
}
