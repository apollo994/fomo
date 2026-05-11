#!/usr/bin/env bash
set -euo pipefail

# Read: <gff_path>\t<count> from the TSV, print only those with count > 0
awk -F'\t' '$2 > 0 {print $1}' ./results/mRNA_count.tsv |
while IFS= read -r gff; do
  dir="$(dirname "$gff")"

  # Find the first .fna in the same directory
  fna="$(find "$dir" -maxdepth 1 -type f -name '*.fna' | head -n 1)"

  # Skip if no fasta found
  [[ -n "${fna:-}" ]] || continue

  # Print the command you want to run
  printf 'bash get_mRNA_seq.sh %q %q\n' "$gff" "$fna"
done
