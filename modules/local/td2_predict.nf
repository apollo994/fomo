process TD2_PREDICT {
    tag "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_low'
    label 'img_td2'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.TD2.pep"), emit: pep
    // TD2 has no --version flag: read the package metadata.
    tuple val("${task.process}"), val('td2'), eval("python -c 'import importlib.metadata as m; print(m.version(\"TD2\"))'"), topic: versions, emit: versions_td2

    script:
    // ext.args  → TD2.LongOrfs flags (default: strand-specific, complete ORFs)
    // ext.args2 → TD2.Predict  flags (default: complete ORFs only)
    def longorfs_args = task.ext.args  ?: '-S --complete-orfs-only'
    def predict_args  = task.ext.args2 ?: '--complete-orfs-only'
    """
    set -euo pipefail

    # Keep the host's ~/.local user site-packages (bind-mounted with \$HOME) off sys.path.
    export PYTHONNOUSERSITE=1

    # TD2.Predict dies when PSAURON has nothing to score, so two cases skip it and write an
    # empty .pep (nothing filtered): an empty input FASTA, and no complete ORF found.
    if [ "\$(grep -c '^>' ${fasta} || true)" -eq 0 ]; then
        : > ${fasta}.TD2.pep
    else
        # `| tee` keeps each step's log and surfaces a traceback in .command.out; pipefail keeps
        # the failure.
        TD2.LongOrfs -t ${fasta} ${longorfs_args} -O td2_work 2>&1 | tee longorfs.log

        n_orfs=\$(grep -c '^>' td2_work/longest_orfs.pep 2>/dev/null || true)
        if [ -z "\${n_orfs}" ] || [ "\${n_orfs}" -eq 0 ]; then
            echo "TD2.LongOrfs found no complete ORFs — skipping TD2.Predict." | tee predict.log
            : > ${fasta}.TD2.pep
        else
            TD2.Predict -t ${fasta} ${predict_args} -O td2_work 2>&1 | tee predict.log
        fi
    fi
    """

    stub:
    """
    touch ${fasta}.TD2.pep
    """
}
