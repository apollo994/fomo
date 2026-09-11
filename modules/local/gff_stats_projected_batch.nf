process GFF_STATS_PROJECTED_BATCH {
    tag "${meta.target_id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_high'

    // Replaces the old per-source GFF_STATS_PROJECTED + GFF_STATS_PROJECTED_TO_MQC pair
    // (S·T tasks each) with ONE task per (target, feature_type, decoy) that loops over
    // every source's projected GFF3 for that target, plus the two aggregate models
    // (allModels_raw, allModels_collapsed). See plans/18_per_target_batching.md — this is
    // the fix for the Nextflow head OOMs the S·T-scaled task count was causing.
    //
    // The JSON gff-feature-stats writes never leaves this task (no longer a separate,
    // published Nextflow output) — gff_stats_to_mqc.py reads it straight back off disk
    // in the same loop iteration.
    //
    // Container is python:3.11 (Debian bookworm, glibc >= 2.34), NOT ubuntu:22.04 as the
    // un-batched GFF_STATS used: this task needs both bin/gff-feature-stats (vendored
    // x86-64 ELF requiring glibc >= 2.34 — see modules/local/gff_stats.nf) and python3
    // (for bin/gff_stats_to_mqc.py). Re-verify the glibc floor if this image tag ever
    // changes: `docker run --rm python:3.11 ldd --version`.
    container 'python:3.11'

    input:
    tuple val(meta), path(gffs)   // every source's *.projected.gff3 for this target, plus
                                   // the *.raw.gff3 and *.collapsed.gff3 aggregate models
                                   // (see subworkflows/local/projection.nf: ch_projected_batch)

    output:
    tuple val(meta), path("*_gffstats_projection_mqc.tsv"), emit: tsv
    tuple val("${task.process}"), val('gff-feature-stats'), eval("gff-feature-stats -V | sed 's/gff-feature-stats //'"), topic: versions, emit: versions_gff_feature_stats
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // Every SPLIT_GFF_BY_SOURCE-produced input file's basename minus '.gff3' is ALREADY
    // the exact string the old per-source ext.prefix built by hand from meta
    // (`<target>.from_<source>.<gtype>.projected` — SPLIT_GFF_BY_SOURCE's --name-template
    // and GFF_STATS_PROJECTED's old ext.prefix already agreed on this, see CLAUDE.md's
    // "merge/split contract"). The two AGGREGATE files do NOT share that convention: they
    // come from FILTER_ALLMODELS/COMBINED_GTF_TO_GFF, published as
    // `<target>.allModels.<gtype>.raw|collapsed.gff3` — the old ext.prefix computed their
    // stats-file name purely from meta (id='allModels_raw'/'allModels_collapsed'),
    // independent of that input filename. So they are symlinked here, under the same
    // '.from_<id>.<gtype>.projected.gff3' convention as the per-source files, using the
    // two known pseudo-source ids — then one glob + one loop handles all S+2 items
    // identically, and no manifest or per-item meta is needed inside the loop itself.
    //
    // The --name given to gff_stats_to_mqc.py is a SEPARATE string from the file prefix:
    // the old ext.sample_name built `<target>.from_<id>.<N>_<class>` (a 1..4 ordering
    // token bin/run_summary_tables.py's own ANNOT_ORDER map parses back out — see there
    // and conf/modules.config's old GFF_STATS_PROJECTED_TO_MQC block). `order` depends only
    // on feature_type/decoy, constant for this whole task, so it is computed once here;
    // only `<id>` varies per file, recovered from its own '.from_<id>.' filename segment.
    def order =
        !meta.decoy && meta.feature_type == 'lnc_RNA' ? '1_lncRNA'        :
         meta.decoy && meta.feature_type == 'lnc_RNA' ? '2_decoy_lncRNA'  :
        !meta.decoy && meta.feature_type == 'mRNA'    ? '3_mRNA'          :
         meta.decoy && meta.feature_type == 'mRNA'    ? '4_decoy_mRNA'    :
                                                        'unknown'
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
        prefix="\${gff%.gff3}"
        # Recover the per-source (or aggregate) id from the '.from_<id>.' segment of the
        # filename itself — this task's single shared meta has no per-file id to read it
        # from. Safe because target/source ids never contain a literal '.' (Ensembl-style
        # names use '_'), so the first '.from_' is unambiguous and stripping to the next
        # '.' isolates exactly the id — the bash equivalent of projection.nf's Groovy
        # string-slice re-key.
        id="\${gff#*.from_}"
        id="\${id%%.*}"
        gff-feature-stats "\$gff" "\${prefix}.stats.json"
        gff_stats_to_mqc.py \\
            --input "\${prefix}.stats.json" \\
            --name "${meta.target_id}.from_\${id}.${order}" \\
            --role projection \\
            --kind projected \\
            --section-id gffstats_projection \\
            --output "\${prefix}_gffstats_projection_mqc.tsv"
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
        prefix="\${gff%.gff3}"
        touch "\${prefix}_gffstats_projection_mqc.tsv"
    done
    """
}
