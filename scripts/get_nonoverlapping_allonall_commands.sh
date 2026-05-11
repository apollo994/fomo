#!/usr/bin/env bash
set -euo pipefail

CMD_DIR="./lncRNA_nonoverlapping_map_commands"
mkdir -p "$CMD_DIR"

# Read lists safely (NUL-delimited)
mapfile -d '' refs < <(find data/ -type f -name '*.fna' -print0)
mapfile -d '' lncs < <(find data/ -type f -name '*.lnc_RNA.longest.nonoverlapping.exons.fa' -print0)

for ref in "${refs[@]}"; do
  # Use the reference parent folder as species name (clean)
  ref_species="$(basename "$(dirname "$ref")")"

  outdir="./results/species/$ref_species"
  mkdir -p "$outdir"

  cmdfile="$CMD_DIR/commands_nonoverlapping_lncmap_${ref_species}.sh"

  for lnc in "${lncs[@]}"; do
    lnc_species="$(basename "$(dirname "$lnc")")"
    bam="$outdir/${ref_species}_lncRNAnonoverlapping_${lnc_species}.bam"

    printf 'bash %q %q %q %q\n' "./run_minimap_base.sh" "$ref" "$lnc" "$bam" >> "$cmdfile"
  done
done
