process GFF_STATS_PROJECTED_BATCH {
    tag "${meta.target_id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_high'

    // One task per (target, feature_type, decoy) over every source's projected GFF3 plus the
    // allModels aggregate (rationale: plans/18). python:3.11 provides both python3 and the
    // glibc >= 2.34 that the vendored gff-feature-stats binary needs.
    container 'python:3.11'

    input:
    tuple val(meta), path(gffs)   // every source's *.projected.gff3 + the *.raw.gff3 aggregate

    output:
    tuple val(meta), path("*_gffstats_projection_mqc.tsv"), emit: tsv
    tuple val("${task.process}"), val('gff-feature-stats'), eval("gff-feature-stats -V | sed 's/gff-feature-stats //'"), topic: versions, emit: versions_gff_feature_stats
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // The aggregate is linked under the per-source '<target>.from_<id>.<track>.projected.gff3'
    // shape so one loop handles every file; the row name is '<target>.from_<id>.<order>', with
    // <id> taken from the filename (ids never contain '.') and <order> fixed per task.
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

    run_one() {
        set -euo pipefail
        gff="\$1"
        prefix="\${gff%.gff3}"
        id="\${gff#*.from_}"
        id="\${id%%.*}"
        gff-feature-stats "\$gff" "\${prefix}.stats.json"
        gff_stats_to_mqc.py \\
            --input "\${prefix}.stats.json" \\
            --name "${meta.target_id}.from_\${id}.${order}" \\
            --role projection \\
            --kind projected \\
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
    for gff in *.projected.gff3; do
        prefix="\${gff%.gff3}"
        touch "\${prefix}_gffstats_projection_mqc.tsv"
    done
    """
}
