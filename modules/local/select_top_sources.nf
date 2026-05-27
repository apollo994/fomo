process SELECT_TOP_SOURCES {
    tag "${meta.id}"
    label 'process_single'

    // Reuse the MultiQC container (has python3); selector is stdlib-only.
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_1'
        : 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'}"

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("*.top_sources.csv"), emit: csv
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    select_top_sources.py \\
        ${stats} \\
        ${args} \\
        --output ${prefix}.top_sources.csv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'source\\n' > ${prefix}.top_sources.csv
    """
}
