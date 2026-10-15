process RENAME_FASTA_HEADERS {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_single'
    label 'img_ubuntu'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.renamed.fasta"), emit: fasta
    // renamed header → original gffread seqname (the GFF3 transcript ID, verbatim): the
    // only place both names are known; FILTER_ANNOTATION matches on it exactly.
    tuple val(meta), path("*.id_map.tsv")   , emit: id_map
    tuple val("${task.process}"), val('mawk'), eval("awk -W version 2>&1 | sed '1!d;s/^mawk //'"), topic: versions, emit: versions_mawk

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    """
    set -euo pipefail

    awk -v type="${meta.label}" -v species="${meta.id}" -v map="${prefix}.id_map.tsv" '
        /^>/ {
            split(\$1, a, ":")
            tid = (length(a) > 1) ? a[2] : substr(\$1, 2)
            print ">" tid "|" type "|" species
            print tid "|" type "|" species "\t" substr(\$1, 2) > map
            next
        }
        { print }
        END { printf "" >> map }   # empty FASTA → empty (but existing) map
    ' ${fasta} > ${prefix}.renamed.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    """
    touch ${prefix}.renamed.fasta ${prefix}.id_map.tsv
    """
}
