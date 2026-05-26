process GFF_TO_GENE_BED {
    tag "${meta.id}"
    label 'process_single'

    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(gff3), path(sizes)

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
        gunzip -c "${gff3}" > _input.gff
    else
        cat "${gff3}" > _input.gff
    fi

    # Sort genes in the order in which chromosomes appear in the genome sizes
    # file (NOT alphabetically) so bedtools complement downstream — which
    # requires BED and genome-file chromosome orders to match — succeeds.
    awk -F'\\t' 'BEGIN{OFS="\\t"} \\
        FNR==NR { order[\$1]=++rank; next } \\
        !/^#/ && NF>=9 && (\$3=="gene"||\$3=="pseudogene"||\$3=="ncRNA_gene") && (\$1 in order) \\
            { print order[\$1], \$1, \$4-1, \$5 }' ${sizes} _input.gff \\
        | sort -k1,1n -k3,3n \\
        | awk -F'\\t' 'BEGIN{OFS="\\t"} \\
            NR==1 { rank=\$1; chrom=\$2; start=\$3; end=\$4; next } \\
            \$1==rank && \$3 <= end { if (\$4 > end) end=\$4; next } \\
            { print chrom, start, end; rank=\$1; chrom=\$2; start=\$3; end=\$4 } \\
            END { if (NR>0) print chrom, start, end }' \\
        > ${prefix}.genes.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genes.bed
    """
}
