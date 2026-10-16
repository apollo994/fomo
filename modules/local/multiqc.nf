process MULTIQC {
    tag "${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"
    label 'process_single'
    label 'img_multiqc'

    // `records` are rendered into custom content first (none for the run summary, whose inputs are
    // custom content already). No versions output: the input waits on the versions topic, so
    // emitting to it would never let the topic close. The report shows its MultiQC version.
    input:
    tuple val(meta), path(records, stageAs: 'records/*'), path(multiqc_files, stageAs: '?/*'), path(multiqc_config, stageAs: '?/*')

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val(meta), path("*_data"), emit: data
    tuple val(meta), path("*_plots"), emit: plots, optional: true

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = [multiqc_config].flatten().collect { cfg -> "--config ${cfg}" }.join(' ')
    def render = records ? "render_mqc.py --records records --target ${meta.target_id} --out ." : ''
    """
    ${render}
    multiqc --force ${args} ${config} ${prefix} .
    """

    stub:
    def prefix = task.ext.prefix ?: 'multiqc_report'
    """
    mkdir ${prefix}_data
    touch ${prefix}_data/.stub ${prefix}.html
    """
}
