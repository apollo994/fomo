#!/usr/bin/env bash
# Part B of ablation/ABLATION_REPORT.md: every configuration side by side, one Nextflow head
# each (all inside this one SLURM job). Each resumes the reference run's cache, so it only
# recomputes what its switched-off component changes. B0 (baseline) is launched separately.
set -uo pipefail
S=/users/rg/fzanarello/pipelines/fomo/ablation/scripts/run_config.sh
$S B1_no_longest       --longest_isoform false                       & sleep 20
$S B2_no_input_td2     -c ablation/configs/B2_no_input_td2.config     & sleep 20
$S B3_no_proj_td2      -c ablation/configs/B3_no_projected_td2.config & sleep 20
$S B4_no_td2           -c ablation/configs/B4_no_td2.config           & sleep 20
$S B5_single_exon      --include_single_exon true                    & sleep 20
$S B6_minimap2_ub      -c ablation/configs/B6_minimap2_ub.config      & sleep 20
$S B7_minimap2_splicehq -c ablation/configs/B7_minimap2_splicehq.config &
wait
