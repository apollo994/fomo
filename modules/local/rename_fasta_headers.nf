process RENAME_FASTA_HEADERS {
    tag "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'

    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.renamed.fasta"), emit: fasta
    // renamed header token → original gffread seqname (= the GFF3 transcript ID=,
    // verbatim). The only place that sees both names, so it is recorded here and
    // FILTER_ANNOTATION builds the filtered GFF3 from it by EXACT id match
    // (plans/20_longest_isoform.md) — no prefix normalisation on either side.
    tuple val(meta), path("*.id_map.tsv")   , emit: id_map
    tuple val("${task.process}"), val('awk'), eval('awk --version 2>&1 | head -n1'), topic: versions, emit: versions_awk

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    def ftype  = meta.feature_type == 'lnc_RNA' ? 'lncRNA' : meta.feature_type
    def type   = meta.decoy ? "decoy_${ftype}" : ftype
    """
    set -euo pipefail

    awk -v type="${type}" -v species="${meta.id}" -v map="${prefix}.id_map.tsv" '
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
