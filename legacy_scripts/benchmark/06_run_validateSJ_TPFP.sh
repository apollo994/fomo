#!/bin/bash
# Run validateSJ on the TP/FP GFF3 splits in benchmark/TPFP_analysis/input/.
# Uses the same Nextflow work dir as scripts/benchmark/02_run_validateSJ.sh so
# STAR alignment results are pulled from the cache via -resume; only the per-
# annotation processes (junction support, transcript support, ...) actually run.
#
# Outputs land in benchmark/TPFP_analysis/validateJS/results/.
set -euo pipefail

ROOT="/users/project/phyloENCODE_009142_no_backup/general_analyses/crossmapping/lycaenidae_benchmark"
WORKDIR="$ROOT/benchmark/TPFP_analysis/validateJS"
PIPELINE="/users/rg/fzanarello/pipelines/validateSJ/main.nf"
GENOME="/users/project/phyloENCODE_009142_no_backup/general_analyses/crossmapping/from_all/data/Polyommatus_icarus_265386_Ensembl_GCA_937595015.1_fe90de4649aacfdb8a5e3d35990743ce/Polyommatus_icarus_265386_Ensembl_GCA_937595015.1_fe90de4649aacfdb8a5e3d35990743ce.1_ilPolIcar1.1_genomic.fna"
NF_WORK="/nfs/scratch01/rg/fzanarello/work_validateSJ"      # SAME as 02_run_validateSJ.sh — required for -resume cache hit
READS_SS="$ROOT/benchmark/validateSJ/reads.csv"             # reuse existing reads samplesheet

cd "$WORKDIR"

nextflow run "$PIPELINE" \
    -profile crg \
    --genome            "$GENOME" \
    --samplesheet       "$WORKDIR/annotations.csv" \
    --samplesheet_reads "$READS_SS" \
    --run_label         lycaenidae_benchmark_TPFP \
    --outdir            "$WORKDIR/results" \
    -work-dir           "$NF_WORK" \
    -resume
