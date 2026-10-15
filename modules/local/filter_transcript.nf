process FILTER_TRANSCRIPT {
    tag "${meta.id}.${meta.feature_type}"
    label 'process_single'
    label 'img_ubuntu'

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*.gff3")                        , emit: gff3
    // One-row funnel (n_transcripts, n_stranded, n_spliced, n_selected), Sample <id>.<ft>.
    tuple val(meta), path("*_transcript_filter_mqc.tsv")   , emit: mqc
    tuple val("${task.process}"), val('mawk'), eval("awk -W version 2>&1 | sed '1!d;s/^mawk //'"), topic: versions, emit: versions_mawk

    script:
    def args = task.ext.args ?: ''            // --longest (params.longest_isoform)
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    // The script reads its input twice, so a gzipped GFF3 is decompressed to node-local
    // scratch ($TMPDIR, then /tmp, then the work dir) and removed on exit.
    """
    if [[ "${gff3}" == *.gz ]]; then
        scratch=\$(mktemp -d "\${TMPDIR:-/tmp}/filter_transcript.XXXXXX" 2>/dev/null \\
            || mktemp -d /tmp/filter_transcript.XXXXXX 2>/dev/null \\
            || mktemp -d "\$PWD/filter_transcript.XXXXXX")
        trap 'rm -rf "\$scratch"' EXIT
        input="\$scratch/input.gff"
        gzip -cd "${gff3}" > "\$input"
    else
        input="${gff3}"
    fi

    filter_transcript.sh \\
        ${args} \\
        --counts ${prefix}_transcript_filter_mqc.tsv \\
        --sample ${meta.id}.${meta.feature_type} \\
        "\$input" ${meta.feature_type} > ${prefix}.gff3
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    """
    touch ${prefix}.gff3
    printf 'Sample\\tfeature\\tn_transcripts\\tn_stranded\\tn_spliced\\tn_selected\\n${meta.id}.${meta.feature_type}\\t${meta.feature_type}\\t0\\t0\\t0\\t0\\n' > ${prefix}_transcript_filter_mqc.tsv
    """
}
