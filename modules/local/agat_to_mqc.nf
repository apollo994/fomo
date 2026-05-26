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
    tuple val(meta), path("*_mqc.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix     ?: "${meta.id}.${meta.kind}"
    def sample     = task.ext.sample_name ?: prefix
    def role       = task.ext.role       ?: 'source'
    def section_id = task.ext.section_id ?: 'agat_input'
    // Filename suffix carries the section so MultiQC's sp patterns can match
    // each section uniquely without overlap (e.g. agat_input vs agat_projection).
    def suffix     = section_id == 'agat_projection' ? '_agat_projection_mqc.tsv'
                                                     : '_agat_input_mqc.tsv'
    """
    agat_to_mqc.py \\
        --input ${yaml} \\
        --name ${sample} \\
        --role ${role} \\
        --section-id ${section_id} \\
        --output ${prefix}${suffix}
    """

    stub:
    def prefix     = task.ext.prefix     ?: "${meta.id}.${meta.kind}"
    def section_id = task.ext.section_id ?: 'agat_input'
    def suffix     = section_id == 'agat_projection' ? '_agat_projection_mqc.tsv'
                                                     : '_agat_input_mqc.tsv'
    """
    touch ${prefix}${suffix}
    """
}
