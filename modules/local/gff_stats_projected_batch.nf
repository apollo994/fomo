process GFF_STATS_PROJECTED_BATCH {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_low'
    label 'img_python'
    label 'amd64_only'

    // One task per (target, track) over every manifest row: each source's projected GFF3 and the
    // allModels aggregate. Each file's identity (its source) comes from the manifest.
    input:
    tuple val(meta), path(gffs), path(manifest)

    output:
    tuple val(meta), path("*_gffstats_projection_mqc.tsv"), emit: tsv
    tuple val(meta), path("*.gffstats.data.json"), emit: data
    tuple val("${task.process}"), val('gff-feature-stats'), eval("gff-feature-stats -V | sed 's/gff-feature-stats //'"), topic: versions, emit: versions_gff_feature_stats
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    // Row name '<target>.from_<source>.<order>_<label>'; the record's identity is the task's
    // plus the row's source.
    """
    set -euo pipefail
    identity='${task.ext.identity}'

    while IFS=\$'\\t' read -r file source; do
        prefix="${meta.target_id}.from_\${source}.${meta.track}.projected"
        gff-feature-stats "\$file" "\${prefix}.stats.json"
        gff_stats_to_mqc.py \\
            --input "\${prefix}.stats.json" \\
            --name "${meta.target_id}.from_\${source}.${meta.order}_${meta.label}" \\
            --role projection \\
            --kind projected \\
            --output "\${prefix}_gffstats_projection_mqc.tsv"
        printf '{"schema":"fomo-data/1","section":"gffstats","identity":{"source":"%s",%s,"data":' "\$source" "\${identity:1}" > "\${prefix}.gffstats.data.json"
        cat "\${prefix}.stats.json" >> "\${prefix}.gffstats.data.json"
        printf '}\\n' >> "\${prefix}.gffstats.data.json"
    done < <(tail -n +2 ${manifest})
    """

    stub:
    """
    tail -n +2 ${manifest} | while IFS=\$'\\t' read -r _file source; do
        prefix="${meta.target_id}.from_\${source}.${meta.track}.projected"
        touch "\${prefix}_gffstats_projection_mqc.tsv" "\${prefix}.gffstats.data.json"
    done
    """
}
