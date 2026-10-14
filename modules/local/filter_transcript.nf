process FILTER_TRANSCRIPT {
    tag "${meta.id}.${meta.feature_type}"
    label 'process_single'

    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(gff3), val(feature_type)

    output:
    tuple val(meta), path("*.gff3")                        , emit: gff3
    // One-row funnel (n_transcripts / n_stranded / n_spliced / n_selected) — MultiQC custom
    // content, and the source of the run summary's dropped_non_longest column.
    // Sample = <id>.<feature_type>: with --include_mrna a species has two rows, and
    // MultiQC keys table rows by sample name, so the bare id would collide.
    tuple val(meta), path("*_transcript_filter_mqc.tsv")   , emit: mqc
    path  "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''            // --longest (params.longest_isoform)
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    // filter_transcript.sh reads its input twice (two awk passes), so a gzipped GFF3
    // needs a real file — decompressed into node-local scratch (SLURM's per-job TMPDIR
    // on the cluster), not the NFS work dir, and removed on exit (plans/26).
    // Fallback order: $TMPDIR, /tmp, the work dir.
    """
    if [[ "${gff3}" == *.gz ]]; then
        scratch=\$(mktemp -d "\${TMPDIR:-/tmp}/filter_transcript.XXXXXX" 2>/dev/null \\
            || mktemp -d /tmp/filter_transcript.XXXXXX 2>/dev/null \\
            || mktemp -d "\$PWD/filter_transcript.XXXXXX")
        trap 'rm -rf "\$scratch"' EXIT
        input="\$scratch/input.gff"
        gzip -cd "${gff3}" > "\$input"
    else
        input="${gff3}"
    fi

    filter_transcript.sh \\
        ${args} \\
        --counts ${prefix}_transcript_filter_mqc.tsv \\
        --sample ${meta.id}.${feature_type} \\
        "\$input" ${feature_type} > ${prefix}.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/[^0-9.]//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}"
    """
    touch ${prefix}.gff3
    printf 'Sample\\tfeature\\tn_transcripts\\tn_stranded\\tn_spliced\\tn_selected\\n${meta.id}.${feature_type}\\t${feature_type}\\t0\\t0\\t0\\t0\\n' > ${prefix}_transcript_filter_mqc.tsv
    touch versions.yml
    """
}
