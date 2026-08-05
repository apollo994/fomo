process GFF_STATS_TO_MQC {
    tag "${meta.id}.${meta.kind}"
    label 'process_single'

    // Stdlib-only JSON parsing — the plain python image is enough (same choice
    // as modules/local/filter_gff_by_id.nf).
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
    def section_id = task.ext.section_id  ?: 'gffstats_input'
    def kind       = meta.kind            ?: 'NA'
    // Each section gets a unique filename suffix so MultiQC's sp patterns match
    // 1:1 with no overlap (gffstats_input vs gffstats_projection vs gffstats_genes).
    def suffix     = section_id == 'gffstats_projection' ? '_gffstats_projection_mqc.tsv'
                                                        : '_gffstats_input_mqc.tsv'
    // The gene-category table is only meaningful for input annotations; projected
    // models carry one synthetic gene per transcript (see BAM_TO_GFF).
    def genes      = task.ext.genes_output ? "--genes-output ${prefix}_gffstats_genes_mqc.tsv" : ''
    """
    gff_stats_to_mqc.py \\
        --input ${json} \\
        --name ${sample} \\
        --role ${role} \\
        --kind ${kind} \\
        --section-id ${section_id} \\
        --output ${prefix}${suffix} \\
        ${genes}
    """

    stub:
    def prefix     = task.ext.prefix     ?: "${meta.id}.${meta.kind}"
    def section_id = task.ext.section_id ?: 'gffstats_input'
    def suffix     = section_id == 'gffstats_projection' ? '_gffstats_projection_mqc.tsv'
                                                        : '_gffstats_input_mqc.tsv'
    def genes      = task.ext.genes_output ? "touch ${prefix}_gffstats_genes_mqc.tsv" : ''
    """
    touch ${prefix}${suffix}
    ${genes}
    """
}
