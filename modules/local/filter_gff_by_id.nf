process FILTER_GFF_BY_ID {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_single'
    label 'img_python'

    // Keep mode (an id map and a FASTA given): keep exactly the transcripts of the FASTA, matched
    // through the id map by exact id; exits 1 if the three disagree. Drop mode (otherwise): drop
    // the listed ids, or nothing when `ids` is []. The GFF3 is staged in a subdirectory so the
    // output may take the same name.
    input:
    tuple val(meta), path(gff3, stageAs: 'input/*'), path(ids), path(fasta)

    output:
    tuple val(meta), path("*.gff3"), emit: gff3
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filter = fasta ? "--keep-fasta ${fasta} --id-map ${ids}" : ids ? "--drop-ids ${ids}" : ''
    """
    filter_gff_by_id.py \\
        --gff ${gff3} \\
        ${filter} \\
        --output ${prefix}.gff3
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff3
    """
}
