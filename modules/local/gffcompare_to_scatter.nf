process GFFCOMPARE_TO_SCATTER {
    tag "${meta.id}"
    label 'process_single'
    label 'img_python'

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("*_mqc.json"), emit: json
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gffcompare_accuracy_mqc.py \\
        ${stats} \\
        --output ${prefix}_mqc.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"id":"gffcompare_accuracy_custom","plot_type":"scatter","data":[]}' > ${prefix}_mqc.json
    """
}
