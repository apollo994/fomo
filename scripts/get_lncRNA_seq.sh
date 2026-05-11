#!/usr/bin/env bash
set -euo pipefail

# ---- args ----
input_gff="${1:?Usage: $0 <annotation.gff3> <genome.fasta>}"
input_fa="${2:?Usage: $0 <annotation.gff3> <genome.fasta>}"

# ---- optional: conda ----
source ~/miniconda3/etc/profile.d/conda.sh
conda activate agat

# ---- filenames (avoid $1"lnc_RNA.gff3" etc.) ----
base="${input_gff%.*}"            # removes last extension (e.g. .gff3)
lnc_gff="${base}.lnc_RNA.gff3"
lnc_longest="${base}.lnc_RNA.longest.gff3"
out_fa="${base}.lnc_RNA.exons.fa"

# ---- tools ----
filter_script="$HOME/dirty_scripts/annotation/filter_transcript.sh"

# ---- run ----
bash "$filter_script" "$input_gff" "lnc_RNA" > "$lnc_gff"
echo "Filter DONE"
agat_sp_keep_longest_isoform.pl -g "$lnc_gff" -o "$lnc_longest"
echo "Longest DONE"
agat_sp_extract_sequences.pl -g "$lnc_longest" -f "$input_fa" -t exon --merge -o "$out_fa"
echo "Extract seq DONE"

echo "Wrote: $out_fa" >&2
