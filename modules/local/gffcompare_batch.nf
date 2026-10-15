process GFFCOMPARE_BATCH {
    tag "${meta.target_id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_high'

    // One task per (target, feature_type, decoy): a separate single-query gffcompare per
    // projected model set (every source, plus the allModels aggregate) against the target's
    // filtered reference. Not gffcompare's multi-query mode, which pools the queries into one
    // Sn/Pr and loses the per-source accuracy. Rationale: plans/18.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffcompare:0.12.6--h9f5acd7_0' :
        'quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0' }"

    input:
    tuple val(meta), path(gffs)          // every source's projected GFF3 + the *.raw.gff3 aggregate
    tuple val(ref_meta), path(reference) // the target's filtered reference GFF3

    output:
    tuple val(meta), path("*.gffcompare.stats"), emit: stats
    tuple val(meta), path("*.tracking"),         emit: tracking
    tuple val(meta), path("*.loci"),             emit: loci
    tuple val(meta), path("*.tmap"),          optional: true, emit: tmap
    tuple val(meta), path("*.refmap"),        optional: true, emit: refmap
    tuple val(meta), path("*.annotated.gtf"), optional: true, emit: annotated_gtf
    tuple val("${task.process}"), val('gffcompare'), eval('gffcompare --version 2>&1 | sed "s/gffcompare v//"'), emit: versions_gffcompare, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Outputs are named after the query file minus '.projected.gff3'; the aggregate is
    // linked under the same '<target>.from_allModels_raw.<track>.projected.gff3' shape.
    // The worker runs in strict mode so a failing gffcompare fails the task.
    """
    shopt -s nullglob

    for gff in *.raw.gff3; do
        ln -s "\$gff" "${meta.target_id}.from_allModels_raw.${meta.feature_type}${meta.decoy ? '.decoy' : ''}.projected.gff3"
    done

    run_one() {
        set -euo pipefail
        gff="\$1"
        prefix="\${gff%.projected.gff3}"
        gffcompare ${args} -r ${reference} -o "\${prefix}.gffcompare" "\$gff"
        mv "\${prefix}.gffcompare" "\${prefix}.gffcompare.stats"
    }
    export -f run_one

    printf '%s\\0' *.projected.gff3 \\
        | xargs -0 -P ${task.cpus} -I{} bash -c 'run_one "\$@"' _ {}
    """

    stub:
    """
    shopt -s nullglob
    for gff in *.raw.gff3; do
        ln -s "\$gff" "${meta.target_id}.from_allModels_raw.${meta.feature_type}${meta.decoy ? '.decoy' : ''}.projected.gff3"
    done
    for gff in *.projected.gff3; do
        prefix="\${gff%.projected.gff3}"
        touch "\${prefix}.gffcompare.stats" "\${prefix}.gffcompare.tracking" "\${prefix}.gffcompare.loci"
    done
    """
}
