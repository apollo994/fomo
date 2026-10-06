process FILTER_ANNOTATION {
    tag "${meta.id}.${meta.feature_type}"
    label 'process_single'

    container 'python:3.11'

    // Builds a source's ONE filtered annotation, <id>.<ft>.filtered.gff3
    // (plans/20_longest_isoform.md): the candidate GFF3 (FILTER_TRANSCRIPT —
    // spliced + longest isoform) subset to EXACTLY the transcripts of the FASTA that
    // is sent to projection (TD2-noncoding for lncRNA, everything for mRNA). The
    // id map from RENAME_FASTA_HEADERS translates each renamed header back to its
    // GFF3 ID verbatim, and filter_gff_by_id.py exits 1 if the three inputs do not
    // describe the same set — so this file is, by construction, what was projected.
    // It is also the gffcompare reference for projections onto this species, and is
    // published under sources/<id>/annotation/.
    //
    // The GFF3 is staged in a subdir (the MAYBE_GUNZIP trick) so ext.prefix can own
    // the full output name without colliding with the input's.
    input:
    tuple val(meta), path(gff3, stageAs: 'input/*'), path(id_map), path(fasta)

    output:
    tuple val(meta), path("*.gff3"), emit: gff3

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}.filtered"
    """
    filter_gff_by_id.py \\
        --gff ${gff3} \\
        --keep-fasta ${fasta} \\
        --id-map ${id_map} \\
        --output ${prefix}.gff3
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}.filtered"
    """
    touch ${prefix}.gff3
    """
}
