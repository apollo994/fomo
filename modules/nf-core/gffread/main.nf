process GFFREAD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), path(gff)
    path fasta

    output:
    tuple val(meta), path("*.gtf")  , emit: gtf             , optional: true
    tuple val(meta), path("*.gff3") , emit: gffread_gff     , optional: true
    tuple val(meta), path("*.fasta"), emit: gffread_fasta   , optional: true
    tuple val(meta), path("*.bed")  , emit: bed             , optional: true
    tuple val("${task.process}"), val('gffread'), eval('gffread --version 2>&1'), topic: versions, emit: versions_gffread

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}"
    def extension   = args.contains("--bed")    ? 'bed' : ( args.contains("-T")       ? 'gtf' : ( ( ['-w', '-x', '-y' ].any { flag -> args.contains(flag) } ) ? 'fasta' : 'gff3' ) )
    def output_name = "${prefix}.${extension}"
    def output      = extension == "fasta"      ? "$output_name" : "-o $output_name"
    def args_sorted = args.replaceAll(/(.*)(-[wxy])(.*)/) { _all, pre, param, post -> "$pre $post $param" }.trim()
    // args_sorted  = Move '-w', '-x', and '-y' to the end of the args string as gffread expects the file name after these parameters
    if ( "$output_name" in [ "$gff", "$fasta" ] ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if (fasta) {
        // FOMO-local deviation from upstream nf-core (plans/26). gffread's `-g` reader has
        // no gzip or bgzf path, and it reads the genome TWICE (a `.fai` scan, then the
        // extraction) — over NFS that was ~2x the uncompressed genome per task, and the
        // reason the decoy/extract tasks stalled on large assemblies. Instead, stage the
        // genome once into node-local scratch (SLURM's per-job TMPDIR on the cluster),
        // decompressing it if gzipped, and run gffread there. The genome therefore stays
        // gzipped from the samplesheet onwards; no GUNZIP task, no plain copy on NFS.
        // `bash -ue` aborts on a failed gzip with its own exit code, so a corrupt archive
        // cannot yield a truncated FASTA. Fallback order: $TMPDIR, /tmp, the work dir.
        def genome = fasta.name.endsWith('.gz') ? fasta.baseName : fasta.name
        def stage  = fasta.name.endsWith('.gz') ? "gzip -cd ${fasta}" : "cat ${fasta}"
        """
        scratch=\$(mktemp -d "\${TMPDIR:-/tmp}/gffread.XXXXXX" 2>/dev/null \\
            || mktemp -d /tmp/gffread.XXXXXX 2>/dev/null \\
            || mktemp -d "\$PWD/gffread.XXXXXX")
        trap 'rm -rf "\$scratch"' EXIT
        ${stage} > "\$scratch/${genome}"

        gffread \\
            $gff \\
            -g "\$scratch/${genome}" \\
            $args_sorted \\
            $output
        """
    } else {
        """
        gffread \\
            $gff \\
            $args_sorted \\
            $output
        """
    }

    stub:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}"
    def extension   = args.contains("--bed")    ? 'bed' : ( args.contains("-T")       ? 'gtf' : ( ( ['-w', '-x', '-y' ].any { flag -> args.contains(flag) } ) ? 'fasta' : 'gff3' ) )
    def output_name = "${prefix}.${extension}"
    if ( "$output_name" in [ "$gff", "$fasta" ] ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch $output_name
    """
}
