process AGAT_TO_MQC {
    tag "${meta.id}.${meta.kind}"
    label 'process_single'

    // Reuse the MultiQC container — it already includes PyYAML (a MultiQC
    // dependency), so we avoid a separate pip install or extra image pull.
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0'
        : 'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(yaml)

    output:
    tuple val(meta), path("*_summary_mqc.tsv"), emit: summary
    tuple val(meta), path("*_full_mqc.tsv"),    emit: full
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    def sample = task.ext.sample_name ?: prefix
    """
    agat_to_mqc.py \\
        --input ${yaml} \\
        --name ${sample} \\
        --summary-output ${prefix}_summary_mqc.tsv \\
        --full-output ${prefix}_full_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    """
    touch ${prefix}_summary_mqc.tsv ${prefix}_full_mqc.tsv
    """
}
