#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
	echo "Usage: $0 <input.bam>" >&2
	exit 1
fi

bam_file="$1"

if [[ ! -f "$bam_file" ]]; then
	echo "Error: input BAM '$bam_file' not found." >&2
	exit 1
fi

tmp_gff="$(mktemp /tmp/convert_bam_to_gff.XXXXXX.gff)"
trap 'rm -f "$tmp_gff"' EXIT

samtools view -h -F 2308 "$bam_file" |                 # exclude unmapped, supplementary and secondary
	bedtools bamtobed -i stdin -bed12 |         # BED12 score column is MAPQ
	awk -v OFS="\t" '{
		# BED12: chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
		split($11, sizes, ",")
		split($12, starts, ",")

		# Output gene and transcript
		print $1, "bedtools", "gene", $2 + 1, $3, $5, $6, ".", "ID=gene-" $4
		print $1, "bedtools", "transcript", $2 + 1, $3, $5, $6, ".", "ID=" $4 ";Parent=gene-" $4


		exon_num = 0
		for (i = 1; i <= int($10); i++) {
				if (sizes[i] <= 0) continue  # This prevents zero size bloks (I between Ns in CIGAR) to create false exons
				exon_num++
				exon_start = $2 + starts[i]
				exon_end = exon_start + sizes[i]
				print $1, "bedtools", "exon", exon_start + 1, exon_end, $5, $6, ".", "ID=" $4 "_exon" exon_num ";Parent=" $4
		}
	}' > "$tmp_gff"


# next block adds NM, AS, de from the bam file into the gff transcript row (9th column)
# and prints the final gff to stdout
printf '##gff-version 3\n'
samtools view -h -F 2308 "$bam_file" |
	awk -v OFS="\t" '
		# First pass: parse SAM and store NM/AS/de by read ID
		NR == FNR {
			if ($0 ~ /^@/) next
			read_id = $1
			if (index(read_id, ";") > 0) {
				print "Error: read name contains ; which would corrupt GFF attributes: " read_id > "/dev/stderr"
				exit 1
			}
			nm = "."
			as = "."
			de = "."
			for (i = 12; i <= NF; i++) {
				if ($i ~ /^NM:i:/) nm = substr($i, 6)
				else if ($i ~ /^AS:i:/) as = substr($i, 6)
				else if ($i ~ /^de:f:/) de = substr($i, 6)
			}
			nm_by_id[read_id] = nm
			as_by_id[read_id] = as
			de_by_id[read_id] = de
			next
		}
		# Second pass: parse GFF and append tags to transcript rows
		{
			if ($3 == "transcript") {
				id = ""
				n = split($9, attrs, ";")
				for (i = 1; i <= n; i++) {
					if (attrs[i] ~ /^ID=/) {
						id = substr(attrs[i], 4)
						break
					}
				}
				nm = (id in nm_by_id ? nm_by_id[id] : ".")
				as = (id in as_by_id ? as_by_id[id] : ".")
				de = (id in de_by_id ? de_by_id[id] : ".")
				$9 = $9 ";NM=" nm ";AS=" as ";de=" de
			}
			print
		}
	' - "$tmp_gff"
