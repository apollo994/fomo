process SAMTOOLS_TO_MQC {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_single'
    label 'img_python'

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("*_samtools_align_mqc.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    def prefix = task.ext.prefix     ?: "${meta.id}"
    def sample = task.ext.sample_name ?: prefix
    """
    samtools_to_mqc.py \\
        --input ${stats} \\
        --name ${sample} \\
        --output ${prefix}_samtools_align_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_samtools_align_mqc.tsv
    """
}
