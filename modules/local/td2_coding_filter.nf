process TD2_CODING_FILTER {
    tag "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'

    container 'python:3.11-slim'

    input:
    tuple val(meta), path(fasta), path(pep)

    output:
    tuple val(meta), path("*.noncoding.fasta"), emit: kept_fasta
    tuple val(meta), path("*.coding_ids.txt"),  emit: coding_ids
    tuple val(meta), path("*_td2_coding_mqc.tsv"), emit: mqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''            // e.g. --psauron-min 0.5
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    def stage  = task.ext.stage  ?: 'input'
    def sample = task.ext.sample_name ?: prefix
    """
    td2_coding_filter.py \\
        --in-fasta ${fasta} \\
        --pep ${pep} \\
        --kept-fasta ${prefix}.noncoding.fasta \\
        --coding-ids ${prefix}.coding_ids.txt \\
        --mqc ${prefix}_td2_coding_mqc.tsv \\
        --sample ${sample} \\
        --stage ${stage} \\
        --feature lncRNA \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    """
    touch ${prefix}.noncoding.fasta ${prefix}.coding_ids.txt ${prefix}_td2_coding_mqc.tsv
    """
}
