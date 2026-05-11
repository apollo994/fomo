#!/usr/bin/env bash
set -euo pipefail

find data -mindepth 1 -maxdepth 1 -type d -print0 |
while IFS= read -r -d '' d; do
  fna="$(find "$d" -maxdepth 1 -type f -name '*.fna' | head -n 1)"
  gff="$(find "$d" -maxdepth 1 -type f \( -name '*.gff' -o -name '*.gff3' -o -name '*.gtf' \) | head -n 1)"

  [[ -n "${fna:-}" && -n "${gff:-}" ]] || continue

  echo "annocli alias" "$gff" "$fna"
done
