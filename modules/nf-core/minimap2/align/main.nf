process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data' :
        'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    val bam_format
    val bam_index_extension
    val cigar_paf_format
    val cigar_bam

    output:
    tuple val(meta), path("*.paf")                       , optional: true, emit: paf
    tuple val(meta), path("*.bam")                       , optional: true, emit: bam
    tuple val(meta), path("*.bam.${bam_index_extension}"), optional: true, emit: index
    tuple val("${task.process}"), val("minimap2"), eval("minimap2 --version"), topic: versions, emit: versions_minimap2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_index = bam_index_extension ? "${prefix}.bam##idx##${prefix}.bam.${bam_index_extension} --write-index" : "${prefix}.bam"
    def cigar_paf = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    def samtools_reset_fastq = bam_input ? "samtools reset --threads ${task.cpus-1} $args3 $reads | samtools fastq --threads ${task.cpus-1} $args4 |" : ''
    def query = bam_input ? "-" : reads
    def target = reference ?: (bam_input ? error("Error: minimap2/align BAM input mode requires reference") : reads)
    if (bam_format && !bam_input) {
        // FOMO-local deviation from upstream nf-core, scoped to the one shape this
        // pipeline actually calls (subworkflows/local/projection.nf: bam_format=true,
        // reads are a FASTA, never sam/bam/cram). Upstream pipes minimap2 straight into
        // `samtools sort` ("-a | samtools sort ..."); without `pipefail`, a killed
        // minimap2 (OOM -> SIGKILL) left samtools sort reading a truncated stream, which
        // reported its OWN generic exit 1 — masking minimap2's real, signal-based exit
        // code (137) that Nextflow's retry errorStrategy keys on
        // (task.exitStatus in (130..145)+104; see conf/base.config). Writing minimap2's
        // output to an intermediate SAM file instead means a killed minimap2 aborts the
        // script (set -e, Nextflow's default .command.sh) with its OWN real exit code
        // before samtools sort ever runs — no pipe, nothing to mask. The bam_input /
        // paf branches below are unused by FOMO and keep the original piped form.
        """
        minimap2 \\
            ${args} \\
            -t ${task.cpus} \\
            ${target} \\
            ${query} \\
            ${cigar_paf} \\
            ${set_cigar_bam} \\
            -a \\
            -o ${prefix}.unsorted.sam

        samtools sort -@ ${task.cpus-1} -o ${bam_index} ${args2} ${prefix}.unsorted.sam
        rm ${prefix}.unsorted.sam
        """
    } else {
        def bam_output = bam_format ? "-a | samtools sort -@ ${task.cpus-1} -o ${bam_index} ${args2}" : "-o ${prefix}.paf"
        """
        $samtools_reset_fastq \\
        minimap2 \\
            ${args} \\
            -t ${task.cpus} \\
            ${target} \\
            ${query} \\
            ${cigar_paf} \\
            ${set_cigar_bam} \\
            ${bam_output}
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_file = bam_format ? "${prefix}.bam" : "${prefix}.paf"
    def bam_index = bam_index_extension ? "touch ${prefix}.bam.${bam_index_extension}" : ""
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    if(bam_input && !reference) {
        error("Error: minimap2/align BAM input mode requires reference!")
	}
    """
    touch $output_file
    ${bam_index}
    """
}
