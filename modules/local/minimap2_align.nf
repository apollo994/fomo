process MINIMAP2_ALIGN {
    tag "${meta.target_id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_medium'
    label 'img_minimap2_samtools'

    input:
    tuple val(meta), path(query), path(index)

    output:
    tuple val(meta), path("*.bam"), path("*.bam.csi"), emit: bam
    tuple val("${task.process}"), val('minimap2'), eval('minimap2 --version'), topic: versions, emit: versions_minimap2
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // meta.exclude_self: the target is also a source, so its own records leave the query and a
    // species is never projected onto itself.
    def reads     = meta.exclude_self ? "${prefix}.query.fa" : query
    def drop_self = meta.exclude_self ? "fasta_drop_species.sh ${meta.target_id} ${query} > ${reads}" : ''
    // minimap2 writes a file rather than a pipe, so a killed minimap2 fails the task with its
    // own exit code (which the retry strategy reads) instead of a truncated sort. CSI, not BAI:
    // BAI cannot address positions above 2^29 bp.
    """
    ${drop_self}

    minimap2 \\
        ${args} \\
        -t ${task.cpus} \\
        ${index} \\
        ${reads} \\
        -a \\
        -o ${prefix}.unsorted.sam

    samtools sort -@ ${task.cpus - 1} -o ${prefix}.bam##idx##${prefix}.bam.csi --write-index ${prefix}.unsorted.sam
    rm -f ${prefix}.unsorted.sam ${drop_self ? reads : ''}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam ${prefix}.bam.csi
    """
}
