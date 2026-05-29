process AGAT_SPSTATISTICS {
    tag "${meta.id}.${meta.kind}"
    label 'process_low'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/agat:1.6.1--pl5321hdfd78af_1'
        : 'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1'}"

    input:
    tuple val(meta), path(gff), path(genome)

    output:
    tuple val(meta), path("*.stats.txt"),  emit: stats_txt
    tuple val(meta), path("*.stats.yaml"), emit: stats_yaml
    tuple val("${task.process}"), val('agat'), eval("agat --version | sed 's/v//'"), topic: versions, emit: versions_agat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    def gs     = genome ? "--gs ${genome}" : ''
    """
    # NOTE: --yaml is a boolean flag in agat_sp_statistics.pl; the YAML
    # filename is derived from --output by appending .yaml. We rename to
    # the cleaner *.stats.yaml form afterwards.
    # --gs <FASTA> enables genome-coverage metrics in the output.
    # if/else (not early exit) so Nextflow's appended eval block always runs.
    if ! awk 'BEGIN{found=0} /^[[:space:]]*#/ {next} NF==0 {next} {found=1; exit} END{exit(found?0:1)}' ${gff}; then
        echo "[WARN] AGAT_SPSTATISTICS: ${gff} has no feature lines; writing empty stats." >&2
        : > ${prefix}.stats.txt
        : > ${prefix}.stats.yaml
    else
        agat_sp_statistics.pl \\
            --gff ${gff} \\
            ${gs} \\
            --output ${prefix}.stats.txt \\
            --yaml \\
            ${args}

        if [[ ! -s ${prefix}.stats.txt.yaml ]]; then
            echo "[ERROR] AGAT_SPSTATISTICS: missing YAML output for ${gff}" >&2
            exit 1
        fi
        mv ${prefix}.stats.txt.yaml ${prefix}.stats.yaml
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    """
    touch ${prefix}.stats.txt ${prefix}.stats.yaml
    """
}
