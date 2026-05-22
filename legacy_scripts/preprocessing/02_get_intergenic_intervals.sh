#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./make_intergenic_bed.sh <annotation.gff3> <out.bed>
#
# Steps:
#  1) Convert ncRNA_gene and pseudogene -> gene
#  2) Add intergenic regions with AGAT
#  3) Extract intergenic_region features and convert to BED6 with awk

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <annotation.gff3> <out.bed>" >&2
  exit 1
fi

# Conda env (as in your original stub)
if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
  # shellcheck disable=SC1090
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
  conda activate agat >/dev/null 2>&1 || true
  unset PERL5LIB PERLLIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
  hash -r
fi

GFF="$1"
OUT_BED="$2"

[[ -s "$GFF" ]] || { echo "ERROR: GFF not found or empty: $GFF" >&2; exit 2; }

command -v agat_sp_add_intergenic_regions.pl >/dev/null || {
  echo "ERROR: agat_sp_add_intergenic_regions.pl not found in PATH (is the 'agat' env active?)." >&2
  exit 3
}

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/agat_intergenic.XXXXXX")"
cleanup() { rm -rf "$tmpdir"; }
trap cleanup EXIT

fixed_gff="${tmpdir}/fixed_types.gff3"
with_intergenic_gff="${tmpdir}/with_intergenic.gff3"

# 1) Convert ncRNA_gene and pseudogene -> gene
awk -F'\t' -v OFS='\t' '
  /^#/ { print; next }
  NF < 9 { print; next }
  {
    if ($3 == "ncRNA_gene" || $3 == "pseudogene") $3 = "gene";
    print
  }
' "$GFF" > "$fixed_gff"

# 2) Add intergenic regions
agat_sp_add_intergenic_regions.pl --gff "$fixed_gff" --out "$with_intergenic_gff"

# 3) Extract intergenic_region and convert to BED6
awk -F'\t' -v OFS='\t' '
  function attr_val(attrs, key,   n,i,a,kv,k,v) {
    n = split(attrs, a, /;/)
    for (i=1; i<=n; i++) {
      split(a[i], kv, /=/)
      k = kv[1]; v = kv[2]
      if (k == key) return v
    }
    return ""
  }
  BEGIN { c=0 }
  /^#/ { next }
  NF < 9 { next }
  $3 != "intergenic_region" { next }
  {
    c++
    chrom=$1
    start0=$4-1
    end=$5
    id=attr_val($9, "ID")
    nm=attr_val($9, "Name")
    name=(id!="" ? id : (nm!="" ? nm : ("intergenic_region_" c)))
    score=0
    strand="."
    print chrom, start0, end, name, score, strand
  }
' "$with_intergenic_gff" > "$OUT_BED"

echo "Wrote: $OUT_BED"
