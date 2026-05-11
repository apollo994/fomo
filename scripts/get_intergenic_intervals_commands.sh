#!/usr/bin/env bash
set -euo pipefail

shopt -s nullglob

for gff in data/*/*.aliasMatch.gff{,3}; do
  out="${gff}.intergenic.bed"
  echo bash get_intergenic_intervals.sh "$gff" "$out"
done
