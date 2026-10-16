process GFF_TO_GENE_BED {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_single'
    label 'img_ubuntu'

    input:
    tuple val(meta), path(gff3), path(sizes)

    output:
    tuple val(meta), path("*.genes.bed"), emit: bed
    tuple val("${task.process}"), val('mawk'), eval("awk -W version 2>&1 | sed '1!d;s/^mawk //'"), topic: versions, emit: versions_mawk

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The GFF3 is streamed (awk's second input); pipefail surfaces a failed gzip.
    def stream = gff3.name.endsWith('.gz') ? "gzip -cd ${gff3}" : "cat ${gff3}"
    """
    set -euo pipefail

    # Genes in the chromosome order of the sizes file, as bedtools complement requires.
    ${stream} \\
        | awk -F'\\t' 'BEGIN{OFS="\\t"} \\
        FNR==NR { order[\$1]=++rank; next } \\
        !/^#/ && NF>=9 && (\$3=="gene"||\$3=="pseudogene"||\$3=="ncRNA_gene") && (\$1 in order) \\
            { print order[\$1], \$1, \$4-1, \$5 }' ${sizes} - \\
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
