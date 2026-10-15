process GFF_STATS_TO_MQC {
    tag "${meta.id}.${meta.kind}"
    label 'process_single'

    container 'python:3.11'

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("*_mqc.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix      ?: "${meta.id}.${meta.kind}"
    def sample     = task.ext.sample_name ?: prefix
    def role       = task.ext.role        ?: 'source'
    def kind       = meta.kind            ?: 'NA'
    def genes      = task.ext.genes_output ? "--genes-output ${prefix}_gffstats_genes_mqc.tsv" : ''
    """
    gff_stats_to_mqc.py \\
        --input ${json} \\
        --name ${sample} \\
        --role ${role} \\
        --kind ${kind} \\
        --output ${prefix}_gffstats_input_mqc.tsv \\
        ${genes}
    """

    stub:
    def prefix     = task.ext.prefix     ?: "${meta.id}.${meta.kind}"
    def genes      = task.ext.genes_output ? "touch ${prefix}_gffstats_genes_mqc.tsv" : ''
    """
    touch ${prefix}_gffstats_input_mqc.tsv
    ${genes}
    """
}
