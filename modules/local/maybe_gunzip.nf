process MAYBE_GUNZIP {
    tag "${archive}"
    label 'process_single'

    // coreutils + gzip (same Wave image the nf-core gunzip module uses)
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(archive, stageAs: 'input/*')   // staged in a subdir so the
                                                          // output can never collide with
                                                          // an uncompressed input of the
                                                          // same name

    output:
    tuple val(meta), path("${output}"), emit: gunzip
    tuple val("${task.process}"), val('gzip'), eval('gzip --version 2>&1 | head -1 | sed "s/^.*) //"'), topic: versions, emit: versions_gzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def isGz          = archive.name.endsWith('.gz')
    def nameWithoutGz = isGz ? archive.baseName : archive.name
    def extension     = file(nameWithoutGz).extension
    def name          = file(nameWithoutGz).baseName
    def prefix        = task.ext.prefix ?: name
    output = prefix + ".${extension}"
    """
    # Decompress if gzipped, otherwise copy through. Either way the output name
    # is derived identically from ext.prefix, so intermediate artefacts are named
    # consistently regardless of whether the input was compressed.
    if [[ "${archive}" == *.gz ]]; then
        gzip -cd ${args} ${archive} > ${output}
    else
        cp -L ${archive} ${output}
    fi
    """

    stub:
    def isGz          = archive.name.endsWith('.gz')
    def nameWithoutGz = isGz ? archive.baseName : archive.name
    def extension     = file(nameWithoutGz).extension
    def name          = file(nameWithoutGz).baseName
    def prefix        = task.ext.prefix ?: name
    output = prefix + ".${extension}"
    """
    touch ${output}
    """
}
