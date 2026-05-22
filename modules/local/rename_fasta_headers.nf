process RENAME_FASTA_HEADERS {
    tag "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'

    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.renamed.fasta"), emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    def ftype  = meta.feature_type == 'lnc_RNA' ? 'lncRNA' : meta.feature_type
    def type   = meta.decoy ? "decoy_${ftype}" : ftype
    """
    awk -v type="${type}" -v species="${meta.id}" '
        /^>/ {
            split(\$1, a, ":")
            tid = (length(a) > 1) ? a[2] : substr(\$1, 2)
            print ">" tid "|" type "|" species
            next
        }
        { print }
    ' ${fasta} > ${prefix}.renamed.fasta
    """
}
