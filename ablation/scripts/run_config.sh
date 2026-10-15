#!/usr/bin/env bash
# Launch one ablation configuration of FOMO (Part B of ablation/ABLATION_REPORT.md).
#
# Every configuration runs in its own launch directory, seeded with a COPY of the cache DB
# of the reference all-vs-all run (session $BASE_SESSION) and resumed from it, against the
# SAME workDir (-profile crg). Only the tasks whose hash the ablation changes are recomputed;
# everything upstream of the switched-off component is reused. The copy keeps the original
# run's launch directory untouched and lets the configurations run side by side.
#
# Usage: run_config.sh NAME [extra nextflow args...]
#   e.g. run_config.sh B1_no_longest --longest_isoform false
#        run_config.sh B4_no_td2 -c ablation/configs/B4_no_td2.config
set -euo pipefail

NAME=$1; shift
REPO=${REPO:-/users/rg/fzanarello/pipelines/fomo}
ABL=${ABL:-/nfs/scratch01/rg/fzanarello/fomo_ablation}
BASE_LAUNCH=${BASE_LAUNCH:-/users/project/phyloENCODE_009142_no_backup/general_analyses/crossmapping/fomo_runs/lycaenidae}
BASE_SESSION=${BASE_SESSION:-3098a02f-6e22-4002-8343-fdf7e8ba9787}
SHEET=${SHEET:-$ABL/samplesheet_both23.csv}

RUN=$ABL/runs/$NAME
mkdir -p "$RUN/.nextflow/cache"
if [[ ! -d "$RUN/.nextflow/cache/$BASE_SESSION" ]]; then
    cp -r "$BASE_LAUNCH/.nextflow/cache/$BASE_SESSION" "$RUN/.nextflow/cache/"
    rm -f "$RUN/.nextflow/cache/$BASE_SESSION/db/LOCK"
fi
# `-resume <id>` looks the session up in .nextflow/history before it opens the cache.
if ! grep -q "$BASE_SESSION" "$RUN/.nextflow/history" 2>/dev/null; then
    grep "$BASE_SESSION" "$BASE_LAUNCH/.nextflow/history" >> "$RUN/.nextflow/history"
fi
cd "$RUN"

# Resolve any relative -c path against the repo, since we cd'd away from it.
args=()
while [[ $# -gt 0 ]]; do
    if [[ $1 == -c ]]; then
        cfg=$2; [[ $cfg = /* ]] || cfg=$REPO/$cfg
        args+=(-c "$cfg"); shift 2
    else
        args+=("$1"); shift
    fi
done

# `-c` is a launcher option, so it has to precede `run`.
cfgs=(); rest=()
for ((i = 0; i < ${#args[@]}; i++)); do
    if [[ ${args[i]} == -c ]]; then cfgs+=(-c "${args[i+1]}"); i=$((i + 1)); else rest+=("${args[i]}"); fi
done

export NXF_OPTS=${NXF_OPTS:--Xms512m -Xmx2g}
nextflow "${cfgs[@]}" run "$REPO" \
    -profile crg \
    --input "$SHEET" \
    --outdir "$RUN/results" \
    --include_decoy \
    -resume "$BASE_SESSION" \
    "${rest[@]}" \
    > "$RUN/nextflow.stdout" 2>&1
