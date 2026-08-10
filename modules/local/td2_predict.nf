process TD2_PREDICT {
    tag "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    // Right-sized from traces: both instances peak at 0.60 GB / 1.24 cores / 15 s against
    // the old process_medium 24 GB / 4 cpu / 1 h (plans/13). PSAURON/TD2 memory scales with
    // transcript count and Lycaenidae lncRNA sets are small, so if a larger clade starts
    // retrying here give TD2_PREDICT its own `withName` block rather than reverting to medium.
    label 'process_low'

    // TD2 (TransDecoder2). BioContainers publishes images only up to td2 1.0.7
    // (no 1.0.8 / 1.1.0 build), so the containers below are Wave builds of the
    // same `bioconda::td2=1.1.0` spec kept in the conda directive for
    // -profile conda. The hashes are derived from that spec — regenerate both
    // (https://seqera.io/containers) if the pin ever changes.
    // PSAURON auto-uses a GPU when torch is CUDA-enabled; these images carry the
    // CPU torch that the spec resolves to, which falls back gracefully. Moving
    // TD2 to a GPU node now means pointing `container` at a CUDA-torch image and
    // adding `--nv` — not swapping the conda spec (see conf/crg.config).
    conda 'bioconda::td2=1.1.0'
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/td2:1.1.0--357f46c1cbbeed20'
        : 'community.wave.seqera.io/library/td2:1.1.0--76046a413a4219c1'}"

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

    # Keep the host's ~/.local/lib/pythonX.Y/site-packages off sys.path. It takes
    # priority over the interpreter's own site-packages whenever the two Python
    # minor versions match, and a stale user-site pandas there shadowed the
    # container/env copy and broke TD2.Predict's import. The container alone does
    # not settle this: singularity/apptainer bind-mounts \$HOME by default
    # (CRG has `mount home = yes`), so user-site is visible inside it too.
    export PYTHONNOUSERSITE=1

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
        # `2>&1 | tee` keeps the per-step log file AND lets the output reach
        # .command.out, so a TD2 traceback shows up in Nextflow's error report
        # instead of being swallowed by the redirect. `set -o pipefail` (above)
        # makes the task still fail when TD2 does, despite tee exiting 0.
        TD2.LongOrfs -t ${fasta} ${longorfs_args} -O td2_work 2>&1 | tee longorfs.log
        TD2.Predict  -t ${fasta} ${predict_args}  -O td2_work 2>&1 | tee predict.log
    fi
    """

    stub:
    """
    touch ${fasta}.TD2.pep
    """
}
