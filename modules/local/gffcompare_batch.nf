process GFFCOMPARE_BATCH {
    tag "${meta.target_id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_high'

    // Replaces the old per-source BENCHMARKING:GFFCOMPARE ((F·S·D+2·F·D)·T_g tasks) with
    // ONE task per (target, feature_type, decoy) that loops S+2 separate single-query
    // gffcompare invocations — every source plus the two aggregate models — against this
    // target's one reference. See plans/18_per_target_batching.md.
    //
    // Deliberately NOT gffcompare's native multi-query mode (`gffcompare -r ref q1 q2 …`):
    // that POOLS every query's transcripts into one shared Sn/Pr, but
    // CONSENSUS_TOP:SELECT_TOP_SOURCES needs a per-source F1 to rank donors
    // (consensus_top.nf) — a pooled comparison would destroy exactly that signal. Each
    // loop iteration is therefore its own complete single-query gffcompare run, under its
    // own `-o` prefix, so its `.stats`/`.tracking`/`.loci`/`.tmap`/`.refmap` are exactly
    // what the un-batched module used to write for that source — only how many times the
    // binary runs per SLURM job has changed.
    conda "${moduleDir}/../nf-core/gffcompare/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffcompare:0.12.6--h9f5acd7_0' :
        'quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0' }"

    input:
    tuple val(meta), path(gffs)          // every source's projected GFF3 for this target,
                                          // plus the *.raw.gff3 / *.collapsed.gff3 aggregates
    tuple val(ref_meta), path(reference) // this target's one filtered reference GFF3
                                          // (BENCHMARKING:FILTER_TARGET) — no fasta (-s)
                                          // input, exactly like the un-batched call, which
                                          // always passed an empty one

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
    // Per-source query files are named `<target>.from_<source>.<gtype>.projected.gff3`
    // (SPLIT_GFF_BY_SOURCE's --name-template), but the OLD un-batched GFFCOMPARE's
    // ext.prefix was `<target>.from_<id>.<gtype>[.decoy].gffcompare` — computed purely
    // from meta, WITHOUT a '.projected' segment, unlike GFF_STATS_PROJECTED's old prefix.
    // So the prefix here strips the longer '.projected.gff3' suffix, not just '.gff3'.
    // The two aggregate files carry no per-file identity of their own at all (see
    // gff_stats_projected_batch.nf's identical note) — symlinked here to the same
    // '.from_<id>.<gtype>.projected.gff3' shape as the per-source files, using the two
    // known pseudo-source ids, so the same suffix-stripping rule applies uniformly.
    """
    shopt -s nullglob

    for gff in *.raw.gff3; do
        ln -s "\$gff" "${meta.target_id}.from_allModels_raw.${meta.feature_type}${meta.decoy ? '.decoy' : ''}.projected.gff3"
    done
    for gff in *.collapsed.gff3; do
        ln -s "\$gff" "${meta.target_id}.from_allModels_collapsed.${meta.feature_type}${meta.decoy ? '.decoy' : ''}.projected.gff3"
    done

    run_one() {
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
    for gff in *.collapsed.gff3; do
        ln -s "\$gff" "${meta.target_id}.from_allModels_collapsed.${meta.feature_type}${meta.decoy ? '.decoy' : ''}.projected.gff3"
    done
    for gff in *.projected.gff3; do
        prefix="\${gff%.projected.gff3}"
        touch "\${prefix}.gffcompare.stats" "\${prefix}.gffcompare.tracking" "\${prefix}.gffcompare.loci"
    done
    """
}
