#!/usr/bin/env bash
set -euo pipefail

# ---- args ----
input_gff="${1:?Usage: $0 <annotation.gff3> <genome.fasta> <feature_type>}"
input_fa="${2:?Usage: $0 <annotation.gff3> <genome.fasta> <feature_type>}"
feature="${3:?Usage: $0 <annotation.gff3> <genome.fasta> <feature_type>}"

# ---- optional: conda ----
source ~/miniconda3/etc/profile.d/conda.sh
conda activate agat

# ---- filenames ----
base="${input_gff%.*}"            # removes last extension (e.g. .gff3)
feat_gff="${base}.${feature}.gff3"
feat_longest="${feat_gff%.*}.longest.gff3"
out_fa="${feat_longest%.*}.exons.fa"

# ---- tools ----
filter_script="$HOME/dirty_scripts/annotation/filter_transcript.sh"

# ---- run ----
bash "$filter_script" "$input_gff" "$feature" > "$feat_gff"
echo "Filter DONE"
agat_sp_keep_longest_isoform.pl -g "$feat_gff" -o "$feat_longest"
echo "Longest DONE"
agat_sp_extract_sequences.pl -g "$feat_longest" -f "$input_fa" -t exon --merge -o "$out_fa"
echo "Extract seq DONE"

echo "Wrote: $out_fa" >&2
