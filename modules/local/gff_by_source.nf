process GFF_BY_SOURCE {
    tag "${meta.target_id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'

    container 'python:3.11'

    // The input GFF3 is staged in a subdir (`stageAs: 'input/*'`) so the `*.gff3`
    // output glob can never pick it up — the outputs are written to the task root and
    // their names are not known until the task runs.
    //
    // One output per source= (SPLIT_GFF_BY_SOURCE in PROJECTION). The subset mode that
    // built the top-N consensus was removed with it (plans/23).
    input:
    tuple val(meta), path(gff3, stageAs: 'input/*')

    // A single source yields a bare Path rather than a List; projection.nf normalises
    // both shapes (`[gffs].flatten()`).
    //
    // optional: true — a target onto which NOTHING projected from ANY
    // source (allModels.raw.gff3 has zero records — a real outcome, not a broken run;
    // seen for real on divergent targets in a large all-vs-all) makes gff_by_source.py
    // write zero files and exit 0 (its own WARNING says so). Without `optional: true`
    // Nextflow raises MissingFileException on the empty glob and fails the task even
    // though the script did exactly what it documents. subworkflows/local/projection.nf
    // handles the resulting empty/missing emit via `join(..., remainder: true)` rather
    // than letting that target silently vanish from the run.
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
