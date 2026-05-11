find . -type f ! -name "*.gz" \( -name "*.aliasMatch.gff" -o -name "*.aliasMatch.gff3" \) -print0 \
| while IFS= read -r -d '' f; do
    n=$(awk -F'\t' '!/^#/ && $3=="mRNA" {c++} END{print c+0}' "$f")
	printf "%s\t%d\n" "$f" "$n"
  done > ./results/mRNA_count.tsv
