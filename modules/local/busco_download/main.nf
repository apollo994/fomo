process BUSCO_DOWNLOAD {
    tag "${lineage}"
    label 'process_single'

    conda "${projectDir}/modules/nf-core/busco/busco/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/6456c1880785adefb4fc9b480bb7662479d5662c17f70d5e8715b7f2a63ee28b/data'
        : 'community.wave.seqera.io/library/busco_numpy:b66937518a305dd7'}"

    input:
    val lineage

    output:
    path "busco_downloads", emit: lineage_path

    script:
    """
    busco --download ${lineage} --download_path busco_downloads
    """
}
