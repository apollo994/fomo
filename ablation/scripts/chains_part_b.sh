#!/usr/bin/env bash
# Chain audit tables for Part B configurations (lncRNA + its decoy track). B0's lncRNA
# outputs are byte-identical to the reference run's, so it reuses offline/full_chains.
# Usage: chains_part_b.sh [CONFIG ...]   (default: every runs/B[1-9]*)
set -euo pipefail
ABL=${ABL:-/nfs/scratch01/rg/fzanarello/fomo_ablation}
S=/users/rg/fzanarello/pipelines/fomo/ablation/scripts/chains_for_run.sh
names=("$@")
[[ ${#names[@]} -gt 0 ]] || names=($(cd "$ABL/runs" && ls -d B[1-9]*))
for n in "${names[@]}"; do
    "$S" "$ABL/runs/$n/results" "$ABL/offline/chains_$n" "$ABL/samplesheet_both23.csv" 'lnc_RNA lnc_RNA.decoy' "${SLURM_CPUS_PER_TASK:-8}"
done
