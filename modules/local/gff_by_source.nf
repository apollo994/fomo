process GFF_BY_SOURCE {
    tag "${meta.target_id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'

    container 'python:3.11'

    // The input GFF3 is staged in a subdir (the MAYBE_GUNZIP trick) so the `*.gff3`
    // output glob can never pick it up — the outputs are written to the task root and
    // their names are not known until the task runs.
    //
    // `keep` is optional: pass [] for split mode (one output per source=), or a
    // one-column `source` CSV for subset mode (a single output holding just those
    // sources). One shape serves both aliases — SPLIT_GFF_BY_SOURCE in PROJECTION and
    // SUBSET_GFF_BY_SOURCE in CONSENSUS_TOP.
    input:
    tuple val(meta), path(gff3, stageAs: 'input/*'), path(keep)

    // Split mode emits N files, subset mode exactly 1 — so this glob yields a List
    // in one case and a bare Path in the other. Both consumers `.transpose()`, which
    // passes a non-List element through unchanged (verified), so the single-file case
    // is not silently dropped and no `arity` declaration is needed.
    output:
    tuple val(meta), path("*.gff3"), emit: gff3

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gff_by_source.py \\
        --gff ${gff3} \\
        --outdir . \\
        ${keep ? "--keep ${keep}" : ''} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff3
    """
}
