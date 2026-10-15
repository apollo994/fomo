process BAM_TO_GFF {
    tag "${meta.target_id}.from_${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'
    label 'img_samtools'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.gff3"), emit: gff3
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam_to_gff.sh ${args} ${bam} > ${prefix}.gff3
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff3
    """
}
