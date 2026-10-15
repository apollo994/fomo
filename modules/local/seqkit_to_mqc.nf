process SEQKIT_TO_MQC {
    tag "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'
    label 'img_ubuntu'

    input:
    tuple val(meta), path(stats_tsv)

    output:
    tuple val(meta), path("*_seqkit_mqc.tsv"), emit: mqc
    tuple val("${task.process}"), val('mawk'), eval("awk -W version 2>&1 | sed '1!d;s/^mawk //'"), topic: versions, emit: versions_mawk

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    def sample = task.ext.sample_name ?: prefix
    """
    set -euo pipefail

    # Rename the first column to Sample and set it to the row name.
    awk -F'\\t' -v OFS='\\t' -v s='${sample}' '
        NR==1 { \$1="Sample"; print; next }
        { \$1=s; print }
    ' ${stats_tsv} > ${prefix}_seqkit_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    """
    touch ${prefix}_seqkit_mqc.tsv
    """
}
