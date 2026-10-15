process TD2_CODING_FILTER {
    tag "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'
    label 'img_python'

    input:
    tuple val(meta), path(fasta), path(pep)

    output:
    tuple val(meta), path("*.noncoding.fasta"), emit: kept_fasta
    tuple val(meta), path("*.coding_ids.txt"),  emit: coding_ids
    tuple val(meta), path("*_td2_coding_mqc.tsv"), emit: mqc
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    def args   = task.ext.args   ?: ''            // e.g. --psauron-min 0.5
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    def stage  = meta.stage ?: 'input'
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
