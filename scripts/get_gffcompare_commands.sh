#!/usr/bin/env bash
set -euo pipefail

RESULTS_DIR="results/species"
DATA_DIR="data"
CMD_DIR="./run_gffcompare_commands"
mkdir -p "$CMD_DIR"

# Read lists safely (NUL-delimited)
mapfile -d '' species < <(find "$RESULTS_DIR" -mindepth 1 -maxdepth 1 -type d -print0)

for sp in "${species[@]}"; do
  echo "WORKING ON" $sp
  ref_species="$(basename "$sp")"
  mapfile -d '' gffs < <(find "$sp" -type f -name "*.gff*" -print0)

  if [[ ${#gffs[@]} -eq 0 ]]; then
    echo "WARNING: no GFF files found for ${ref_species}" >&2
    continue
  fi

  cmdfile="$CMD_DIR/run_gffcompare_commands_${ref_species}.sh"
  : > "$cmdfile"

  for gff in "${gffs[@]}"; do

	if [[ "$gff" == *"_lncRNAnonoverlapping_"* ]]; then
		ref_type="aliasMatch.lnc_RNA.longest.nonoverlapping.gff3"
	elif [[ "$gff" == *"_lncRNA_"* ]]; then
		ref_type="aliasMatch.lnc_RNA.longest.gff3"
	elif [[ "$gff" == *"_decoy_"* ]]; then
		ref_type="aliasMatch.lnc_RNA.longest.nonoverlapping.decoy.gff3"
	elif [[ "$gff" == *"_mRNA_"* ]]; then
		ref_type="aliasMatch.mRNA.gff3"
	else
		echo "WARNING: unrecognized gene type for ${gff}, skipping" >&2
		continue
	fi

	ref_file=$(find "${DATA_DIR}/${ref_species}" -maxdepth 1 -name "${ref_species}*.${ref_type}" 2>/dev/null | head -n1 || true)

	if [[ ! -f "$ref_file" ]]; then
		echo "WARNING: ref file not found for ${ref_species} (${ref_type}), skipping" >&2
		continue
	fi

	out_dir="$(dirname "$gff")/$(basename "$gff" .gff).gffcompare"
	printf 'mkdir -p %q && gffcompare %q -r %q -o %q/gffcompare\n' \
		"$out_dir" "$gff" "$ref_file" "$out_dir" >> "$cmdfile"
  done
done
