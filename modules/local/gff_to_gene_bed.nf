process GFF_TO_GENE_BED {
    tag "${meta.id}"
    label 'process_single'

    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*.genes.bed"), emit: bed

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ "${gff3}" == *.gz ]]; then
        gunzip -c "${gff3}"
    else
        cat "${gff3}"
    fi | awk -F'\\t' 'BEGIN{OFS="\\t"} !/^#/ && NF>=9 && (\$3=="gene"||\$3=="pseudogene"||\$3=="ncRNA_gene"){print \$1,\$4-1,\$5}' \\
       | sort -k1,1 -k2,2n \\
       > ${prefix}.genes.bed
    """
}
