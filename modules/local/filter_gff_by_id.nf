process FILTER_GFF_BY_ID {
    tag "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'

    container 'python:3.11'

    input:
    tuple val(meta), path(gff3), path(drop_ids)

    output:
    tuple val(meta), path("*.filtered.gff3"), emit: gff3

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    """
    filter_gff_by_id.py \\
        --gff ${gff3} \\
        --drop-ids ${drop_ids} \\
        --output ${prefix}.filtered.gff3
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    """
    touch ${prefix}.filtered.gff3
    """
}
