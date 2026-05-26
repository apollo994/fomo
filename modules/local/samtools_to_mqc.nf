process SAMTOOLS_TO_MQC {
    tag "${meta.id}"
    label 'process_single'

    // Pure Python — reuse the MultiQC container (Python + standard libs).
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_1'
        : 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'}"

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("*_samtools_align_mqc.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

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
