process GFF_STATS {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_low'
    label 'img_ubuntu'
    label 'amd64_only'

    // One task per GFF3, or per (target, track) batch of projected GFF3s with their manifest
    // (file, source). Each record wraps the gff-feature-stats JSON verbatim; its identity is the
    // task's, plus the row's source. bin/gff-feature-stats (provenance:
    // bin/gff-feature-stats.PROVENANCE) reads .gff3.gz; its input must be grouped by seqid.
    input:
    tuple val(meta), path(gffs), path(manifest)

    output:
    tuple val(meta), path("*.gffstats.data.json"), emit: data
    tuple val("${task.process}"), val('gff-feature-stats'), eval("gff-feature-stats -V | sed 's/gff-feature-stats //'"), topic: versions, emit: versions_gff_feature_stats

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rows   = manifest ? "tail -n +2 ${manifest}" : "printf '%s\\t\\n' ${gffs}"
    """
    set -euo pipefail
    identity='${task.ext.identity}'

    while IFS=\$'\\t' read -r file source; do
        name="${prefix}"
        record_identity="\$identity"
        if [[ -n "\$source" ]]; then
            name="${meta.target_id}.from_\${source}.${meta.track}.projected"
            record_identity="{\\"source\\":\\"\${source}\\",\${identity:1}"
        fi
        gff-feature-stats "\$file" "\${name}.stats.json"
        printf '{"schema":"fomo-data/1","section":"gffstats","identity":%s,"data":' "\$record_identity" > "\${name}.gffstats.data.json"
        cat "\${name}.stats.json" >> "\${name}.gffstats.data.json"
        printf '}\\n' >> "\${name}.gffstats.data.json"
    done < <(${rows})
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rows   = manifest ? "tail -n +2 ${manifest}" : "printf '%s\\t\\n' ${gffs}"
    """
    ${rows} | while IFS=\$'\\t' read -r _file source; do
        name="${prefix}"
        [[ -n "\$source" ]] && name="${meta.target_id}.from_\${source}.${meta.track}.projected"
        echo '{}' > "\${name}.gffstats.data.json"
    done
    """
}
