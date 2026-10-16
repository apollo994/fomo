process MINIMAP2_INDEX {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_index'
    label 'img_minimap2'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi"), emit: index
    tuple val("${task.process}"), val('minimap2'), eval('minimap2 --version'), topic: versions, emit: versions_minimap2

    script:
    def args = task.ext.args ?: ''
    """
    minimap2 \\
        -t ${task.cpus} \\
        -d ${fasta.baseName}.mmi \\
        ${args} \\
        ${fasta}
    """

    stub:
    """
    touch ${fasta.baseName}.mmi
    """
}
