process BEDTOOLS_COMPLEMENT {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_single'
    label 'img_bedtools'

    input:
    tuple val(meta), path(bed), path(sizes)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed 's/bedtools v//'"), topic: versions, emit: versions_bedtools

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.complement"
    """
    bedtools complement -i ${bed} -g ${sizes} ${args} > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.complement"
    """
    touch ${prefix}.bed
    """
}
