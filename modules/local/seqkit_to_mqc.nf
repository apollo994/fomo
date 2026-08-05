process SEQKIT_TO_MQC {
    tag "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'

    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(stats_tsv)

    output:
    tuple val(meta), path("*_seqkit_mqc.tsv"), emit: mqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    def sample = task.ext.sample_name ?: prefix
    """
    set -euo pipefail

    # Rewrite the seqkit TSV so MultiQC's custom_content table uses a clean
    # sample identifier (matching GFF stats row names) instead of the long FASTA
    # filename. The first column is renamed from 'file' to 'Sample' and the
    # value is replaced with \$sample on every data row.
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
