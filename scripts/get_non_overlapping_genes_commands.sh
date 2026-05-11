#!/usr/bin/env bash
set -euo pipefail

for spdir in ./data/*/; do
  sp="$(basename "$spdir")"

  # Match files (adjust patterns if needed)
  lnc_files=("$spdir"/*lnc_RNA.longest.gff3)
  pcg_files=("$spdir"/*mRNA.longest.gff3)

  # Skip if not found
  [[ -e "${lnc_files[0]}" ]] || { echo "[WARN] No lncRNA.longest.gff3 for $sp" >&2; continue; }
  [[ -e "${pcg_files[0]}" ]] || { echo "[WARN] No mRNA.longest.gff3 for $sp" >&2; continue; }

  lnc="${lnc_files[0]}"
  pcg="${pcg_files[0]}"

  base="$(basename "$lnc" .gff3)"
  out="$spdir/${base}.nonoverlapping.gff3"

  echo bash ./get_non_overlapping_genes.sh "$lnc" "$pcg" "$out"
done
