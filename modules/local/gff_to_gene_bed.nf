process GFF_TO_GENE_BED {
    tag "${meta.id}"
    label 'process_single'

    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*.genes.bed"), emit: bed
    tuple val("${task.process}"), val('awk'), eval('awk --version 2>&1 | head -n1'), topic: versions, emit: versions_awk

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail

    if [[ "${gff3}" == *.gz ]]; then
        gunzip -c "${gff3}"
    else
        cat "${gff3}"
    fi | awk -F'\\t' 'BEGIN{OFS="\\t"} !/^#/ && NF>=9 && (\$3=="gene"||\$3=="pseudogene"||\$3=="ncRNA_gene"){print \$1,\$4-1,\$5}' \\
       | sort -k1,1 -k2,2n \\
       | awk -F'\\t' 'BEGIN{OFS="\\t"} \\
           NR==1 { chrom=\$1; start=\$2; end=\$3; next } \\
           \$1==chrom && \$2 <= end { if (\$3 > end) end=\$3; next } \\
           { print chrom, start, end; chrom=\$1; start=\$2; end=\$3 } \\
           END { if (NR>0) print chrom, start, end }' \\
       > ${prefix}.genes.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genes.bed
    """
}
