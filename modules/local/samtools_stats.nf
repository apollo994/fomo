process SAMTOOLS_STATS {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_single'
    label 'img_samtools'

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools stats ${args} --threads ${task.cpus} ${bam} > ${prefix}.stats
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats
    """
}
