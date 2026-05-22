#!/bin/bash
# Run the validateSJ Nextflow pipeline against all GFF3 annotations in
# benchmark/input/ using paired-end RNA-seq from RNAseq/.
#
# Results land in benchmark/validateSJ/results/.
# Nextflow work dir lives on scratch (not backed up).
set -euo pipefail

ROOT="/users/project/phyloENCODE_009142_no_backup/general_analyses/crossmapping/lycaenidae_benchmark"
WORKDIR="$ROOT/benchmark/validateSJ"
PIPELINE="/users/rg/fzanarello/pipelines/validateSJ/main.nf"
GENOME="/users/project/phyloENCODE_009142_no_backup/general_analyses/crossmapping/from_all/data/Polyommatus_icarus_265386_Ensembl_GCA_937595015.1_fe90de4649aacfdb8a5e3d35990743ce/Polyommatus_icarus_265386_Ensembl_GCA_937595015.1_fe90de4649aacfdb8a5e3d35990743ce.1_ilPolIcar1.1_genomic.fna"
NF_WORK="/nfs/scratch01/rg/fzanarello/work_validateSJ"

cd "$WORKDIR"

nextflow run "$PIPELINE" \
    -profile crg \
    --genome            "$GENOME" \
    --samplesheet       "$WORKDIR/annotations.csv" \
    --samplesheet_reads "$WORKDIR/reads.csv" \
    --run_label         lycaenidae_benchmark \
    --outdir            "$WORKDIR/results" \
    -work-dir           "$NF_WORK" \
    -resume
