#!/usr/bin/env bash
set -euo pipefail

RESULTS_DIR="results/species"
DATA_DIR="data"
CMD_DIR="./run_gffread_commands"
mkdir -p "$CMD_DIR"

# Read species dirs safely (NUL-delimited)
mapfile -d '' species < <(find "$RESULTS_DIR" -mindepth 1 -maxdepth 1 -type d -print0)

for sp in "${species[@]}"; do
  ref_species="$(basename "$sp")"

  # Each subdirectory is a query species whose genome was used for crossmapping
  mapfile -d '' query_dirs < <(find "$sp" -mindepth 1 -maxdepth 1 -type d -print0)

  if [[ ${#query_dirs[@]} -eq 0 ]]; then
    echo "WARNING: no query species subdirs found under ${ref_species}" >&2
    continue
  fi

  # Find the ref species genome FASTA once, shared across all query dirs
  genome=$(find "${DATA_DIR}/${ref_species}" -maxdepth 1 -name "*.fna" 2>/dev/null | head -n1 || true)
  if [[ ! -f "$genome" ]]; then
    echo "WARNING: genome FASTA not found for ${ref_species}, skipping" >&2
    continue
  fi

  cmdfile="$CMD_DIR/run_gffread_commands_${ref_species}.sh"
  : > "$cmdfile"

  for qdir in "${query_dirs[@]}"; do

    mapfile -d '' gffs < <(find "$qdir" -maxdepth 1 -type f -name "*.gff" -print0)
    if [[ ${#gffs[@]} -eq 0 ]]; then
      echo "WARNING: no GFF files found in ${qdir}, skipping" >&2
      continue
    fi

    for gff in "${gffs[@]}"; do
      out_fa="$(dirname "$gff")/$(basename "$gff" .gff).mrna.fa"
      printf 'gffread %q -g %q -w %q\n' "$gff" "$genome" "$out_fa" >> "$cmdfile"
    done
  done
done
