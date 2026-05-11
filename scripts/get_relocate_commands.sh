#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_relocation_all.sh /path/to/relocate.py
#
# Expects, for each species folder under data/*/ :
#   - one intergenic BED:  *.intergenic.bed
#   - one nonoverlapping GFF/GFF3: *.nonoverlapping.gff or *.nonoverlapping.gff3
#
# Writes, in the same species folder:
#   <nonoverlapping_basename>.relocated.lncRNA.gff3

PY="./relocate_loci2.py"

command -v python3 >/dev/null || { echo "ERROR: python3 not found" >&2; exit 2; }
[[ -f "$PY" ]] || { echo "ERROR: relocate script not found: $PY" >&2; exit 2; }

shopt -s nullglob

for dir in data/*/; do
  beds=( "$dir"/*.intergenic.bed )
  gffs=( "$dir"/*.aliasMatch.lnc_RNA.longest.nonoverlapping.gff3 )

  # require exactly one of each (skip otherwise)
  if (( ${#beds[@]} != 1 || ${#gffs[@]} != 1 )); then
    echo "Skipping $dir (found ${#beds[@]} intergenic bed(s), ${#gffs[@]} nonoverlapping gff(s))" >&2
    continue
  fi

  bed="${beds[0]}"
  gff="${gffs[0]}"

  base="$(basename "$gff")"
  base_noext="${base%.*}"   # removes .gff or .gff3
  out="${dir}/${base_noext}.decoy.gff3"

  # echo "[$(basename "$dir")]"
  # echo "  BED: $bed"
  # echo "  GFF: $gff"
  # echo "  OUT: $out"

  echo python3 "$PY" "$bed" "$gff" "$out"

done
