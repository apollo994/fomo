process GFFCOMPARE_TO_SCATTER {
    tag "${meta.id}"
    label 'process_single'

    // Reuse the MultiQC container (has python3); parser is stdlib-only.
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_1'
        : 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'}"

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("*_mqc.json"), emit: json
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

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
