process EXTRACT_TRANSCRIPTS {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_low'
    label 'img_gffread'

    input:
    tuple val(meta), path(gff), path(genome)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val("${task.process}"), val('gffread'), eval('gffread --version 2>&1'), topic: versions, emit: versions_gffread

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def plain  = genome.name.endsWith('.gz') ? genome.baseName : genome.name
    def stage  = genome.name.endsWith('.gz') ? "gzip -cd ${genome}" : "cat ${genome}"
    // gffread reads the genome twice and cannot read gzip, so the genome is staged once, plain,
    // on node-local scratch ($TMPDIR, then /tmp, then the work dir) and removed on exit.
    """
    scratch=\$(mktemp -d "\${TMPDIR:-/tmp}/gffread.XXXXXX" 2>/dev/null \\
        || mktemp -d /tmp/gffread.XXXXXX 2>/dev/null \\
        || mktemp -d "\$PWD/gffread.XXXXXX")
    trap 'rm -rf "\$scratch"' EXIT
    ${stage} > "\$scratch/${plain}"

    gffread ${gff} -g "\$scratch/${plain}" -w ${prefix}.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    """
}
