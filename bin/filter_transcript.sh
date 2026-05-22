#!/usr/bin/env bash
# filter_transcript.sh
# Usage: bash filter_transcript.sh input.gff[3] TRANSCRIPT_TYPE > output.gff
set -euo pipefail

gff="${1:-}"
ttype="${2:-}"

if [[ -z "${gff}" || -z "${ttype}" ]]; then
  echo "Usage: bash $0 input.gff[3] TRANSCRIPT_TYPE > output.gff" >&2
  exit 1
fi

awk -F'\t' -v TTYPE="$ttype" '
  BEGIN { OFS="\t" }

  # ---------- PASS 1: collect transcript IDs (of type TTYPE) and their Parent IDs ----------
  FNR==NR {
    if ($0 ~ /^#/ || NF < 9) next
    if ($3 == TTYPE) {
      id=""; parent=""
      n=split($9,a,";")
      for (i=1;i<=n;i++) {
        if (a[i] ~ /^ID=/)      { id=a[i]; sub(/^ID=/,"",id) }
        else if (a[i] ~ /^Parent=/) { parent=a[i]; sub(/^Parent=/,"",parent) }
      }
      if (id != "") tx[id]=1

      # Parent can be comma-separated (or absent)
      if (parent != "") {
        m=split(parent,p,",")
        for (j=1;j<=m;j++) if (p[j]!="") parent_of_tx[p[j]]=1
      }
    }
    next
  }

  # ---------- helpers ----------
  function get_attr(field, key,   n,i,a,kv) {
    n=split(field,a,";")
    for(i=1;i<=n;i++){
      if (a[i] ~ ("^" key "=")) {
        kv=a[i]
        sub("^" key "=", "", kv)
        return kv
      }
    }
    return ""
  }

  # ---------- PASS 2: print lines matching your rules ----------
  /^#/ { next }
  NF < 9 { next }

  {
    id     = get_attr($9, "ID")
    parent = get_attr($9, "Parent")

    # Rule 1: print if this feature is a parent of a selected transcript
    # (works for gene, nc_gene, or any other parent feature type)
    if (id != "" && parent_of_tx[id]) { print; next }

    # Rule 2: print if this feature ID is one of the selected transcripts (the ttype lines)
    if (id != "" && tx[id]) { print; next }

    # Rule 3: print if this feature is a child of a selected transcript (exon, CDS, UTR, etc.)
    if (parent != "") {
      m=split(parent,p,",")
      for (j=1;j<=m;j++) {
        if (tx[p[j]]) { print; next }
      }
    }
  }
' "$gff" "$gff"
