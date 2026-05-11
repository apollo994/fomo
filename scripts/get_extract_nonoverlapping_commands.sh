#!/usr/bin/env bash
set -euo pipefail



find ./data/ -type f -name '*.lnc_RNA.longest.nonoverlapping.gff3' -print0 |
while IFS= read -r -d '' gff; do
  base="${gff%.*}"                 # remove last extension (.gff3)
  out_fa="${base}.exons.fa"

  dir="$(dirname "$gff")"

  # Find the first .fna in the same directory (if multiple, this picks the first)
  fna="$(find "$dir" -maxdepth 1 -type f -name '*.fna' -print | head -n 1)"

  if [[ -z "$fna" ]]; then
    echo "[WARN] No .fna found next to: $gff" >&2
    continue
  fi
  echo "source ~/miniconda3/etc/profile.d/conda.sh && conda activate agat && agat_sp_extract_sequences.pl -g \"$gff\" -f \"$fna\" -t exon --merge -o \"$out_fa\""
done
