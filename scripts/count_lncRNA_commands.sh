find . -type f ! -name "*.gz" \( -name "*.aliasMatch.gff" -o -name "*.aliasMatch.gff3" \) -print0 \
| while IFS= read -r -d '' f; do
    n=$(awk -F'\t' '!/^#/ && $3=="lnc_RNA" {c++} END{print c+0}' "$f")
	printf "%s\t%d\n" "$f" "$n"
  done > ./results/lnc_RNA_count.tsv
