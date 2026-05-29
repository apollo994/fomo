process TD2_PREDICT {
    tag "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/td2:1.0.7--pyhdfd78af_0'
        : 'quay.io/biocontainers/td2:1.0.7--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.TD2.pep"), emit: pep
    tuple val("${task.process}"), val('td2'), eval('TD2.Predict --version 2>&1 | sed "s/.* //"'), topic: versions, emit: versions_td2

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args  → TD2.LongOrfs flags (default: strand-specific, complete ORFs)
    // ext.args2 → TD2.Predict  flags (default: complete ORFs only)
    def longorfs_args = task.ext.args  ?: '-S --complete-orfs-only'
    def predict_args  = task.ext.args2 ?: '--complete-orfs-only'
    """
    set -euo pipefail

    # LongOrfs intermediates go in td2_work/; TD2.Predict reads that via -O and
    # writes <fasta>.TD2.pep to the current working directory (the task root),
    # so the *.TD2.pep output glob matches directly.
    TD2.LongOrfs -t ${fasta} ${longorfs_args} -O td2_work > longorfs.log 2>&1
    TD2.Predict  -t ${fasta} ${predict_args}  -O td2_work > predict.log  2>&1
    """

    stub:
    """
    touch ${fasta}.TD2.pep
    """
}
