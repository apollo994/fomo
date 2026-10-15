process MERGE_FASTA {
    tag "${meta.id}"
    label 'process_single'

    // coreutils only
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    // Concatenate the per-source spliced FASTAs of one gene type into the all-sources FASTA
    // that minimap2 aligns in one run. Records carry their species (`>tid|type|species`), so
    // the merge is lossless. The caller sorts the input list: the task hash keys on its order.
    input:
    tuple val(meta), path(fastas, stageAs: 'input/*')

    output:
    tuple val(meta), path("*.fasta"), emit: fasta

    when:
    task.ext.when == null || task.ext.when

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
