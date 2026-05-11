#!/usr/bin/env bash
set -euo pipefail

# ---- args ----
input_gff="${1:?Usage: $0 <annotation.gff3> <genome.fasta>}"
input_fa="${2:?Usage: $0 <annotation.gff3> <genome.fasta>}"

# ---- optional: conda ----
source ~/miniconda3/etc/profile.d/conda.sh
conda activate agat

# ---- filenames ----
base="${input_gff%.*}"            # removes last extension (e.g. .gff3)
mRNA_gff="${base}.mRNA.gff3"
mRNA_longest="${base}.mRNA.longest.gff3"
out_fa="${base}.mRNA.exons.fa"

# ---- tools ----
filter_script="$HOME/dirty_scripts/annotation/filter_transcript.sh"

# ---- run ----
bash "$filter_script" "$input_gff" "mRNA" > "$mRNA_gff"
echo "Filter DONE"
agat_sp_keep_longest_isoform.pl -g "$mRNA_gff" -o "$mRNA_longest"
echo "Longest DONE"
agat_sp_extract_sequences.pl -g "$mRNA_longest" -f "$input_fa" -t exon --merge -o "$out_fa"
echo "Extract seq DONE"

echo "Wrote: $out_fa" >&2
