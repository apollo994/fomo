find ./data -type f ! -name "*.gz" \( -name "*.aliasMatch.lnc_RNA.longest.nonoverlapping.gff3" -o -name "*.aliasMatch.lnc_RNA.longest.nonoverlapping.gff" \) -print0 \
| while IFS= read -r -d '' f; do
    n=$(awk -F'\t' '!/^#/ && $3=="lnc_RNA" {c++} END{print c+0}' "$f")
	printf "%s\t%d\n" "$f" "$n"
  done > ./results/lnc_RNA_nonoverlapping_count.tsv
