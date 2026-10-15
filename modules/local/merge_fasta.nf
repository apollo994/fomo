process MERGE_FASTA {
    tag "${meta.id}"
    label 'process_single'
    label 'img_ubuntu'

    // Concatenate the per-source spliced FASTAs of one gene type into the all-sources FASTA
    // that minimap2 aligns in one run. Records carry their species (`>tid|type|species`), so
    // the merge is lossless. The caller sorts the input list: the task hash keys on its order.
    input:
    tuple val(meta), path(fastas, stageAs: 'input/*')

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val("${task.process}"), val('coreutils'), eval("cat --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_coreutils

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail

    cat ${fastas} > ${prefix}.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    """
}
