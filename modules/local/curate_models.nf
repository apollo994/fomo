process CURATE_MODELS {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_low'
    label 'img_curate'

    // `ref` is the target's samplesheet GFF3 (every biotype), or [] for a target without one.
    input:
    tuple val(meta), path(gff3, stageAs: 'input/*'), path(ref, stageAs: 'ref/*')

    // curated:  <target>.curated.<ft>[.decoy].gff3.gz         — always (may hold zero genes)
    // merged:   <target>.curated.<ft>[.decoy].merged.gff3.gz  — reference + curated, only with a reference
    // report:   <prefix>.curation.{tsv,json}                  — observed / discarded / novel
    // mqc:      three *_curation_{observed,discarded,novel}_mqc.tsv (one row each)
    output:
    tuple val(meta), path("${prefix}.gff3.gz")                 , emit: curated
    tuple val(meta), path("${prefix}.merged.gff3.gz")          , emit: merged, optional: true
    tuple val(meta), path("${prefix}.curation.{tsv,json}")     , emit: report
    tuple val(meta), path("${prefix}_curation_*_mqc.tsv")      , emit: mqc
    tuple val(meta), path("${prefix}.curation.data.json")      , emit: data
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('gffcompare'), eval('gffcompare --version 2>&1 | sed "s/gffcompare v//"'), topic: versions, emit: versions_gffcompare
    tuple val("${task.process}"), val('gffread'), eval('gffread --version 2>&1'), topic: versions, emit: versions_gffread

    script:
    // No `def`: the output declarations read `prefix`.
    prefix = task.ext.prefix ?: "${meta.target_id}.curated.${meta.track}"
    def args = task.ext.args ?: ''
    """
    curate_models.py \\
        --gff ${gff3} \\
        ${ref ? "--ref ${ref}" : ''} \\
        --prefix ${prefix} \\
        --target ${meta.target_id} \\
        --track ${meta.track} \\
        --identity '${task.ext.identity}' \\
        ${args}

    curation_to_mqc.py \\
        --json ${prefix}.curation.json \\
        --prefix ${prefix} \\
        --sample ${meta.target_id}.${meta.track}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.target_id}.curated.${meta.track}"
    """
    python3 -c "import gzip, sys; [gzip.open(p, 'wt').close() for p in sys.argv[1:]]" ${prefix}.gff3.gz ${ref ? "${prefix}.merged.gff3.gz" : ''}
    touch ${prefix}.curation.tsv ${prefix}.curation.json ${prefix}.curation.data.json
    for t in observed discarded novel; do touch ${prefix}_curation_\${t}_mqc.tsv; done
    """
}
