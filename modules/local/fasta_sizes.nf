process FASTA_SIZES {
    tag "${meta.id}"
    label 'process_single'

    // Chromosome sizes for the decoy chain, in one streaming pass over the (still gzipped)
    // genome. Same output as `cut -f1,2` of a .fai; samtools faidx is not used because it
    // rejects plain gzip.
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.sizes"), emit: sizes
    tuple val("${task.process}"), val('awk'), eval('awk --version 2>&1 | head -n1'), topic: versions, emit: versions_awk

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stream = fasta.name.endsWith('.gz') ? "gzip -cd ${fasta}" : "cat ${fasta}"
    """
    set -euo pipefail

    ${stream} \\
        | LC_ALL=C awk 'BEGIN { OFS = "\\t" }
            /^>/ { if (n != "") print n, len; n = substr(\$1, 2); len = 0; next }
                 { sub(/\\r\$/, ""); len += length(\$0) }
            END  { if (n != "") print n, len }' \\
        > ${prefix}.sizes
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sizes
    """
}
