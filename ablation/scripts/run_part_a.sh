#!/usr/bin/env bash
# Part A of ablation/ABLATION_REPORT.md: offline ablations on the reference all-vs-all run.
# Needs chains_for_run.sh output in $ABL/offline/full_chains first.
set -euo pipefail
ABL=${ABL:-/nfs/scratch01/rg/fzanarello/fomo_ablation}
FULL=${FULL:-/users/project/phyloENCODE_009142_no_backup/general_analyses/crossmapping/fomo_runs/lycaenidae/results}
PY="python3 /users/rg/fzanarello/pipelines/fomo/ablation/scripts/ablate.py"
SHEET=$ABL/samplesheet_both23.csv
cd "$ABL/offline"
$PY sweep     --sheet "$SHEET" --results "$FULL" --chains full_chains --out A1_sweep.tsv --feature-types lnc_RNA mRNA --kmax 12 &
$PY calib     --sheet "$SHEET" --results "$FULL" --chains full_chains --out A3_calib.tsv --feature-types lnc_RNA mRNA &
$PY donors    --sheet "$SHEET" --results "$FULL" --chains full_chains --out A4_donors.tsv --reps 5 &
$PY persource --sheet "$SHEET" --results "$FULL" --out A5_persource.tsv &
wait
