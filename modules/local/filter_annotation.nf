process FILTER_ANNOTATION {
    tag "${meta.id}.${meta.feature_type}"
    label 'process_single'

    container 'python:3.11'

    // A source's filtered annotation: the candidate GFF3 subset to exactly the transcripts of
    // the projected FASTA, matched through the id map by exact id; exits 1 if the three
    // disagree. Rationale: plans/20. The GFF3 is staged in a subdirectory so the output may
    // take the same name.
    input:
    tuple val(meta), path(gff3, stageAs: 'input/*'), path(id_map), path(fasta)

    output:
    tuple val(meta), path("*.gff3"), emit: gff3

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}.filtered"
    """
    filter_gff_by_id.py \\
        --gff ${gff3} \\
        --keep-fasta ${fasta} \\
        --id-map ${id_map} \\
        --output ${prefix}.gff3
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}.filtered"
    """
    touch ${prefix}.gff3
    """
}
