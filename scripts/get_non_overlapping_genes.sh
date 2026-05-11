#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <gff1.gff3> <gff2.gff3> <out.gff3>" >&2
  exit 1
fi

module load BEDTools
source ~/miniconda3/etc/profile.d/conda.sh
conda activate agat

GFF1="$1"
GFF2="$2"
OUT="$3"

# Exit early if output already exists and is non-empty
if [[ -s "$OUT" ]]; then
  echo "[INFO] Output already exists (non-empty): $OUT"
  exit 0
fi

# Minimal tool checks
command -v bedtools >/dev/null || { echo "[ERROR] bedtools not in PATH" >&2; exit 2; }
command -v agat_sp_filter_feature_from_keep_list.pl >/dev/null || { echo "[ERROR] agat_sp_filter_feature_from_keep_list.pl not in PATH" >&2; exit 2; }

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

g1_gene_bed="$tmpdir/gff1.genes.bed"
g2_gene_bed="$tmpdir/gff2.features.bed"
g1_gene_bed_s="$tmpdir/gff1.genes.sorted.bed"
g2_gene_bed_s="$tmpdir/gff2.features.sorted.bed"
nonoverlap_genes_bed="$tmpdir/gff1.genes.nonoverlap_any_gff2_feature.bed"
keep_ids="$tmpdir/genes_to_keep.ids.txt"

# echo "[INFO] Keeping genes (and children) from: $(basename "$GFF1") that do NOT overlap ANY feature in: $(basename "$GFF2")"

# 1) Extract gene-like features from GFF1 -> BED (gene/pseudogene/ncRNA_gene)
awk -F'\t' 'BEGIN{OFS="\t"}
  $0 !~ /^#/ && ($3=="gene" || $3=="pseudogene" || $3=="ncRNA_gene") {
    id=""
    if (match($9, /ID=([^;]+)/, m)) id=m[1]
    if (id=="") next
    print $1, $4-1, $5, id, 0, $7
  }' "$GFF1" > "$g1_gene_bed"

# 2) Extract ALL features from GFF2 -> BED (any feature type counts as overlap)
awk -F'\t' 'BEGIN{OFS="\t"}
  $0 !~ /^#/ && ($3=="gene" || $3=="pseudogene" || $3=="ncRNA_gene") {
    id=""
    if (match($9, /ID=([^;]+)/, m)) id=m[1]
    if (id=="") next
    print $1, $4-1, $5, id, 0, $7
  }' "$GFF2" > "$g2_gene_bed"


echo "[INFO] gff1 ${GFF1}: $(wc -l < "$g1_gene_bed") genes"
echo "[INFO] gff2 ${GFF2}: $(wc -l < "$g2_gene_bed") genes"


# 3) Sort BED
sort -k1,1 -k2,2n "$g1_gene_bed" > "$g1_gene_bed_s"
sort -k1,1 -k2,2n "$g2_gene_bed" > "$g2_gene_bed_s"

# 4) Keep only genes from GFF1 with ZERO overlap against ANY feature in GFF2
# Add -s if you want strand-specific "no overlap"
bedtools intersect -v \
  -a "$g1_gene_bed_s" \
  -b "$g2_gene_bed_s" \
  > "$nonoverlap_genes_bed" || true

# 5) If none remain, write empty output (valid GFF header) and exit
if [[ ! -s "$nonoverlap_genes_bed" ]]; then
  printf "##gff-version 3\n" > "$OUT"
  echo "[INFO] No non-overlapping genes found. Wrote empty GFF3: $OUT"
  exit 0
fi

# 6) Build keep list of gene IDs
cut -f4 "$nonoverlap_genes_bed" | sort -u > "$keep_ids"

echo "[INFO] Keep list: $(wc -l < "$keep_ids") genes"

# 7) Keep only those genes (and their children) from GFF1
agat_sp_filter_feature_from_keep_list.pl \
  --gff "$GFF1" \
  --keep_list "$keep_ids" \
  -o "$OUT" >/dev/null

out_genes_n=$(awk -F'\t' '$0 !~ /^#/ && ($3=="gene" || $3=="pseudogene" || $3=="ncRNA_gene"){c++} END{print c+0}' "$OUT")
removed_genes_n=$(( $(wc -l < "$g1_gene_bed") - out_genes_n ))

echo -e "\n[INFO] Done. Wrote: $OUT \n[INFO] output_genes=${out_genes_n}  removed_genes=${removed_genes_n}"

# Final check
if [[ ! -s "$OUT" ]]; then
  echo "[ERROR] Output not created or empty: $OUT" >&2
  exit 5
fi

# Remove AGAT report derived from output name
base="${OUT%.*}"
agat_report="${base}_report.txt"
rm -f -- "$agat_report"




