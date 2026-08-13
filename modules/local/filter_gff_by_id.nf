process FILTER_GFF_BY_ID {
    tag "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'

    container 'python:3.11'

    // The GFF3 is staged in a subdir (the MAYBE_GUNZIP trick) so ext.prefix can own
    // the full output name — including the case where output and input would
    // otherwise be called the same thing.
    input:
    tuple val(meta), path(gff3, stageAs: 'input/*'), path(drop_ids)

    output:
    tuple val(meta), path("*.gff3"), emit: gff3

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    """
    filter_gff_by_id.py \\
        --gff ${gff3} \\
        --drop-ids ${drop_ids} \\
        --output ${prefix}.gff3
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    """
    touch ${prefix}.gff3
    """
}
