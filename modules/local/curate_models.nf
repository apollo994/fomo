process CURATE_MODELS {
    tag "${meta.target_id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    // Phase-A real-data runs (plans/21): median 4 s, max 28 s, peak RSS 772 MB on the largest of
    // 82 targets (191 k models from 176 species, 1.14 M reference exons) — 1 core, well inside 4 GB.
    label 'process_low'

    // bin/curate_models.py needs python3 (stdlib only) + gffcompare (exact intron-chain cliques)
    // + gffread (validates the curated GFF3) in ONE task. No BioContainers image has all three,
    // so these are Seqera Containers (Wave) builds of the conda spec below — the gffcompare /
    // gffread pins match modules/nf-core/{gffcompare,gffread} so curation and benchmarking run
    // the same versions. Regenerate both (https://seqera.io/containers) if a pin changes.
    conda 'conda-forge::python=3.11 bioconda::gffcompare=0.12.6 bioconda::gffread=0.12.7'
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gffcompare_gffread:970e511c83295d93'
        : 'community.wave.seqera.io/library/python_gffcompare_gffread:8f0a261560fe9c3b'}"

    // `ref` is the target's samplesheet GFF3 — the FULL annotation, every biotype — or [] for a
    // pure target. Both inputs are staged in subdirs so the output globs can never pick them up.
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
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('gffcompare'), eval('gffcompare --version 2>&1 | sed "s/gffcompare v//"'), topic: versions, emit: versions_gffcompare
    tuple val("${task.process}"), val('gffread'), eval('gffread --version 2>&1'), topic: versions, emit: versions_gffread

    when:
    task.ext.when == null || task.ext.when

    script:
    // No `def`: the output declarations above read `prefix` (same idiom as nf-core gunzip).
    prefix = task.ext.prefix ?: "${meta.target_id}.curated.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    def args = task.ext.args ?: ''
    def track = "${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    """
    curate_models.py \\
        --gff ${gff3} \\
        ${ref ? "--ref ${ref}" : ''} \\
        --prefix ${prefix} \\
        --target ${meta.target_id} \\
        --track ${track} \\
        ${args}

    curation_to_mqc.py \\
        --json ${prefix}.curation.json \\
        --prefix ${prefix} \\
        --sample ${meta.target_id}.${track}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.target_id}.curated.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    """
    echo | gzip > ${prefix}.gff3.gz
    ${ref ? "echo | gzip > ${prefix}.merged.gff3.gz" : ''}
    touch ${prefix}.curation.tsv ${prefix}.curation.json
    for t in observed discarded novel; do touch ${prefix}_curation_\${t}_mqc.tsv; done
    """
}
