process RELOCATE_LOCI {
    tag "${meta.id}.${meta.feature_type}"
    label 'process_single'

    container 'python:3.11-slim'

    input:
    tuple val(meta), path(gff3), path(intergenic_bed)

    output:
    tuple val(meta), path("*.decoy.gff3"), emit: gff3

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    def seed   = (params.relocate_seed != null) ? "--seed ${params.relocate_seed}" : ''
    """
    relocate_loci.py \\
        --input-gff ${gff3} \\
        --intergenic-bed ${intergenic_bed} \\
        --output-gff ${prefix}.decoy.gff3 \\
        --feature-type ${meta.feature_type} \\
        ${seed} \\
        ${args}
    """
}
