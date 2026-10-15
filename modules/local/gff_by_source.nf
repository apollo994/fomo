process GFF_BY_SOURCE {
    tag "${meta.target_id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'

    container 'python:3.11'

    // One output per source= value; the input is staged in a subdirectory so the output
    // glob cannot match it.
    input:
    tuple val(meta), path(gff3, stageAs: 'input/*')

    // A single source yields a bare Path, not a List. An input with zero records (nothing
    // projected onto the target) writes no file: the output is optional.
    output:
    tuple val(meta), path("*.gff3"), optional: true, emit: gff3

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gff_by_source.py \\
        --gff ${gff3} \\
        --outdir . \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def template = args =~ /--name-template\s+'([^']+)'/
    def name = template.find() ? template.group(1).replace('{source}', 'stub') : "${meta.id}.gff3"
    """
    touch ${name}
    """
}
