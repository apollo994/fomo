process TD2_PREDICT {
    tag "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_medium'

    // TD2 (TransDecoder2) runs via conda (no container build needed). td2 is a
    // noarch conda package, so it resolves on amd64/arm64/osx-arm64 alike.
    // PSAURON auto-uses a GPU when the env's torch is CUDA-enabled — the CRG
    // profile overrides this spec with conda-forge::pytorch-gpu + a GPU SLURM
    // request; the default CPU torch falls back gracefully everywhere else.
    conda 'bioconda::td2=1.1.0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.TD2.pep"), emit: pep
    // TD2 has no --version flag; read the installed version from package metadata
    // (verified in the bioconda td2=1.1.0 conda env).
    tuple val("${task.process}"), val('td2'), eval("python -c 'import importlib.metadata as m; print(m.version(\"TD2\"))'"), topic: versions, emit: versions_td2

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args  → TD2.LongOrfs flags (default: strand-specific, complete ORFs)
    // ext.args2 → TD2.Predict  flags (default: complete ORFs only)
    def longorfs_args = task.ext.args  ?: '-S --complete-orfs-only'
    def predict_args  = task.ext.args2 ?: '--complete-orfs-only'
    """
    set -euo pipefail

    # Guard against an empty input FASTA (e.g. a distant source whose lncRNA do
    # not project onto the target → 0 transcripts). With no sequences TD2.LongOrfs
    # writes no psauron_score.csv and TD2.Predict crashes reading it, so emit an
    # empty .pep and skip — downstream the coding-ID set is then simply empty.
    if [ "\$(grep -c '^>' ${fasta} || true)" -eq 0 ]; then
        : > ${fasta}.TD2.pep
    else
        # LongOrfs intermediates go in td2_work/; TD2.Predict reads that via -O and
        # writes <fasta>.TD2.pep to the current working directory (the task root),
        # so the *.TD2.pep output glob matches directly.
        TD2.LongOrfs -t ${fasta} ${longorfs_args} -O td2_work > longorfs.log 2>&1
        TD2.Predict  -t ${fasta} ${predict_args}  -O td2_work > predict.log  2>&1
    fi
    """

    stub:
    """
    touch ${fasta}.TD2.pep
    """
}
