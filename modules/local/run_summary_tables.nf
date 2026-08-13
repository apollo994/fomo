process RUN_SUMMARY_TABLES {
    tag 'run_summary'
    label 'process_single'

    // Stdlib-only TSV/JSON parsing — the plain python image is enough (same choice
    // as modules/local/gff_stats_to_mqc.nf). Deliberately NOT the MultiQC image:
    // that container is pinned in one place (modules/nf-core/multiqc) and this
    // process needs nothing from it.
    container 'python:3.11'

    input:
    // The pipeline-wide union of MultiQC inputs, staged into one flat directory.
    // Every basename is unique today (target-scoped prefixes all start with
    // ${meta.target_id}., source-side ones are species-keyed), so a collision here
    // is Nextflow telling you a naming contract broke — which is the point.
    path(mqc_files, stageAs: 'mqc/*')
    // May be EMPTY: no annotated target ⇒ no ranking ⇒ no top_sources.csv.
    path(top_sources, stageAs: 'top_sources/*')
    path roles_csv

    output:
    path("*_mqc.{tsv,json}"), emit: mqc
    path("run_summary.json"), emit: tables
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // `?:` on the List guards the empty-channel case; the instanceof guard keeps a
    // single-file input (a bare Path, not a List) working.
    def top_list = top_sources instanceof List ? top_sources : (top_sources ? [top_sources] : [])
    def tops = top_list ? "--top-sources ${top_list.join(' ')}" : ''
    """
    run_summary_tables.py \\
        --mqc-dir mqc \\
        --roles ${roles_csv} \\
        ${tops} \\
        ${args}
    """

    stub:
    """
    printf 'Metric\\tValue\\nn_species\\t0\\n' > run_overview_mqc.tsv
    echo '{}' > run_summary.json
    """
}
