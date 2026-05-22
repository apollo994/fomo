process FILTER_TRANSCRIPT {
    tag "${meta.id}.${meta.feature_type}"
    label 'process_single'

    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(gff3), val(feature_type)

    output:
    tuple val(meta), path("*.gff3"), emit: gff3
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    """
    if [[ "${gff3}" == *.gz ]]; then
        gunzip -c "${gff3}" > _input.gff
    else
        ln -sf "${gff3}" _input.gff
    fi

    filter_transcript.sh _input.gff ${feature_type} > ${prefix}.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/[^0-9.]//g')
    END_VERSIONS
    """
}
