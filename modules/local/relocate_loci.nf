process RELOCATE_LOCI {
    tag "${meta.id}.${meta.feature_type}"
    label 'process_single'

    container 'python:3.11-slim'

    input:
    tuple val(meta), path(gff3), path(intergenic_bed)

    output:
    tuple val(meta), path("*.decoy.gff3"), emit: gff3
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    def seed       = (params.relocate_seed != null) ? "--seed ${params.relocate_seed}" : ''
    def decoy_cap  = (params.decoy_cap     != null) ? "--decoy-cap ${params.decoy_cap}" : ''
    """
    relocate_loci.py \\
        --input-gff ${gff3} \\
        --intergenic-bed ${intergenic_bed} \\
        --output-gff ${prefix}.decoy.gff3 \\
        --feature-type ${meta.feature_type} \\
        ${seed} \\
        ${decoy_cap} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    """
    touch ${prefix}.decoy.gff3
    """
}
