#!/usr/bin/env bash
set -euo pipefail

CMD_DIR="./mRNA_map_commands"
mkdir -p "$CMD_DIR"

# Read lists safely (NUL-delimited)
mapfile -d '' refs < <(find data/ -type f -name '*.fna' -print0)
mapfile -d '' mRNAs < <(find data/ -type f -name '*mRNA.exons.fa' -print0)

for ref in "${refs[@]}"; do
  # Use the reference parent folder as species name (clean)
  ref_species="$(basename "$(dirname "$ref")")"

  outdir="./results/species/$ref_species"
  mkdir -p "$outdir"

  cmdfile="$CMD_DIR/commands_mRNAmap_${ref_species}.sh"

  for mRNA in "${mRNAs[@]}"; do
    mRNA_species="$(basename "$(dirname "$mRNA")")"
    bam="$outdir/${ref_species}_mRNA_${mRNA_species}.bam"

    printf 'bash %q %q %q %q\n' "./run_minimap_base.sh" "$ref" "$mRNA" "$bam" >> "$cmdfile"
  done
done
