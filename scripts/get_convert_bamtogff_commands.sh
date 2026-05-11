#!/usr/bin/env bash
set -euo pipefail

CMD_DIR="./convert_bamtogff_commands"
mkdir -p "$CMD_DIR"

# Read lists safely (NUL-delimited)
mapfile -d '' species < <(find ./results/species -mindepth 1 -maxdepth 1 -type d -print0)

for sp in "${species[@]}"; do
  ref_species="$(basename "$sp")"
  mapfile -d '' bams < <(find "$sp" -type f -name "*.bam" -print0)
  echo $sp
  cmdfile="$CMD_DIR/convert_bamtogff_commands_${ref_species}.sh"
  : > "$cmdfile"

  for bam in "${bams[@]}"; do
    gff="$(dirname "$bam")/$(basename "$bam" .bam).gff"
    printf 'bash %q %q > %q\n' "../scripts/convert_bam_to_gff.sh" "$bam" "$gff" >> "$cmdfile"
  done
done
