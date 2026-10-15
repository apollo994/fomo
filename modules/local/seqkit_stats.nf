process SEQKIT_STATS {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_low'
    label 'img_seqkit'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: stats
    tuple val(meta), path("*.fasta_stats.data.json"), emit: data
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), topic: versions, emit: versions_seqkit

    script:
    def args   = task.ext.args ?: '--all'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqkit stats --tabular --threads ${task.cpus} ${args} ${fasta} > ${prefix}.tsv
    seqkit_stats_record.sh '${task.ext.identity}' ${prefix}.tsv > ${prefix}.fasta_stats.data.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv ${prefix}.fasta_stats.data.json
    """
}
