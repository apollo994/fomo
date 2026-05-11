#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

for sp_dir in ./data/*/; do
  sp="${sp_dir%/}"
  sp="${sp##*/}"

  # collect gff/gff3 files
  mapfile -t gffs < <(find "./data/$sp" -type f \( -name "*.nonoverlapping.gff3" -o -name "*.nonoverlapping.gff" \) -print)

  # pick a reference .fna (first match). Adjust if you expect multiple.
  ref="$(find "./data/$sp" -type f -name "*.fna" -print -quit || true)"

  if [[ -z "${ref}" ]]; then
    echo "WARN: no .fna found for species '$sp' in ./data/$sp (skipping)" >&2
    continue
  fi

  if [[ "${#gffs[@]}" -eq 0 ]]; then
    echo "WARN: no .gff/.gff3 found for species '$sp' in ./data/$sp (skipping)" >&2
    continue
  fi

  for gff in "${gffs[@]}"; do
    gff_base="$(basename "$gff")"
    out="./data/$sp/${gff_base}.stats"

    echo "Running: $gff  (ref: $ref) -> $out" >&2
    bash ~/dirty_scripts/annotation/extract_features.sh "$gff" "$ref" > "$out"
  done
done
