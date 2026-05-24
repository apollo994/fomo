process AGAT_SPSTATISTICS {
    tag "${meta.id}.${meta.kind}"
    label 'process_low'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/agat:1.6.1--pl5321hdfd78af_1'
        : 'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1'}"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.stats.txt"),  emit: stats_txt
    tuple val(meta), path("*.stats.yaml"), emit: stats_yaml
    tuple val("${task.process}"), val('agat'), eval("agat --version | sed 's/v//'"), topic: versions, emit: versions_agat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    """
    # NOTE: --yaml is a boolean flag in agat_sp_statistics.pl; the YAML
    # filename is derived from --output by appending .yaml. We rename to
    # the cleaner *.stats.yaml form afterwards.
    agat_sp_statistics.pl \\
        --gff ${gff} \\
        --output ${prefix}.stats.txt \\
        --yaml \\
        ${args}

    mv ${prefix}.stats.txt.yaml ${prefix}.stats.yaml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    """
    touch ${prefix}.stats.txt ${prefix}.stats.yaml
    """
}
