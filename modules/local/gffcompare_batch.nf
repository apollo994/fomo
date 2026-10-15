process GFFCOMPARE_BATCH {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_low'
    label 'img_gffcompare'

    // One task per (target, track): one gffcompare run per manifest row (every source and the
    // allModels aggregate) against the target's filtered reference. Rationale: plans/18.
    input:
    tuple val(meta), path(gffs), path(manifest), path(reference)

    output:
    tuple val(meta), path("*.gffcompare.stats"), emit: stats
    tuple val(meta), path("*.tracking"),         emit: tracking
    tuple val(meta), path("*.loci"),             emit: loci
    tuple val(meta), path("*.tmap"),          optional: true, emit: tmap
    tuple val(meta), path("*.refmap"),        optional: true, emit: refmap
    tuple val(meta), path("*.annotated.gtf"), optional: true, emit: annotated_gtf
    tuple val(meta), path("*.accuracy.data.json"), emit: data
    tuple val("${task.process}"), val('gffcompare'), eval('gffcompare --version 2>&1 | sed "s/gffcompare v//"'), emit: versions_gffcompare, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    gffcompare_batch.sh \\
        --manifest ${manifest} \\
        --reference ${reference} \\
        --target ${meta.target_id} \\
        --track ${meta.track} \\
        --identity '${task.ext.identity}' \\
        -- ${args}
    """

    stub:
    """
    tail -n +2 ${manifest} | while IFS=\$'\\t' read -r _file source; do
        prefix="${meta.target_id}.from_\${source}.${meta.track}"
        touch "\${prefix}.gffcompare.stats" "\${prefix}.gffcompare.tracking" "\${prefix}.gffcompare.loci" "\${prefix}.accuracy.data.json"
    done
    """
}
