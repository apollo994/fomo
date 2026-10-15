process RELOCATE_LOCI {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_single'
    label 'img_python'

    input:
    tuple val(meta), path(gff3), path(intergenic_bed)

    output:
    tuple val(meta), path("*.decoy.gff3"), emit: gff3
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    """
    relocate_loci.py \\
        --input-gff ${gff3} \\
        --intergenic-bed ${intergenic_bed} \\
        --output-gff ${prefix}.decoy.gff3 \\
        --feature-type ${meta.feature_type} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    """
    touch ${prefix}.decoy.gff3
    """
}
