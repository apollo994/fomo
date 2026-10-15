process RUN_SUMMARY_TABLES {
    tag 'run_summary'
    label 'process_single'

    container 'python:3.11'

    input:
    // The pipeline-wide union of MultiQC inputs in one flat directory: basenames are unique
    // (target- or species-prefixed), so a collision means a naming contract broke.
    path(mqc_files, stageAs: 'mqc/*')
    path roles_csv

    output:
    path("*_mqc.{tsv,json}"), emit: mqc
    path("run_summary.json"), emit: tables
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    run_summary_tables.py \\
        --mqc-dir mqc \\
        --roles ${roles_csv} \\
        ${args}
    """

    stub:
    """
    printf 'Metric\\tValue\\nn_species\\t0\\n' > run_overview_mqc.tsv
    echo '{}' > run_summary.json
    """
}
