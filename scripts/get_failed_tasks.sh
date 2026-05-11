#!/usr/bin/env bash
set -euo pipefail

LOG_GLOB="/nfs/scratch01/rg/fzanarello/logs/allmRNA_21064690_*.out"
PATTERN="All steps completed"
THRESH=123
RUN_CMDS="./run_allonall_mRNA.sh"

[[ -f "$RUN_CMDS" ]] || { echo "ERROR: cannot find $RUN_CMDS" >&2; exit 1; }

grep -c -- "$PATTERN" $LOG_GLOB 2>/dev/null \
  | awk -F: -v thr="$THRESH" '$2 < thr {print $1 ":" $2}' \
  | while IFS=: read -r logfile count; do
      task_id="$(basename "$logfile" | sed -E 's/.*_([0-9]+)\.out/\1/')"
      [[ "$task_id" =~ ^[0-9]+$ ]] || continue

      # zero-based task_id -> one-based line number for sed
      line_no=$((task_id + 1))

      cmd="$(sed -n "${line_no}p" "$RUN_CMDS")"
      [[ -n "$cmd" ]] || continue

      printf "%s\t%s\t%s\n" "$task_id" "$count" "$cmd"
    done
