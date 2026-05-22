#!/bin/bash
#
# Generates minimap2 alignment commands for every (target genome, source transcript-set)
# combination across all species in ../../input/. For each ordered pair (target != source),
# three jobs are emitted: mRNA -> pcg/, lnc_RNA (nonoverlapping) -> lnc/, decoy -> decoy/.
#
# It also (re)creates the result folder hierarchy under ../../minimap_models/, mirroring the
# subdirectory layout used in ../../liftoff_models/ (template: <target>_from_<source>/{pcg,lnc,decoy}).
#
# Usage:
#   bash ./00_run_minimap_generator.sh                        # print commands to stdout
#   bash ./00_run_minimap_generator.sh > minimap_tasks.txt    # capture for slurm array
#   bash ./submit.sh minimap_tasks.txt minimap                # launch via send_array.sh

set -euo pipefail

INPUT_DIR="../../input"
MINIMAP_DIR="../../minimap_models"

TEMPLATE=$'pcg\nlnc\ndecoy'

species=( $(ls "$INPUT_DIR") )

# Build target_from_source/{pcg,lnc,decoy} hierarchy and emit commands
for target in "${species[@]}"; do
	target_fna=$(find "$INPUT_DIR/$target" -name "*.fna")
	for source in "${species[@]}"; do
		[[ "$target" == "$source" ]] && continue
		pair="${target}_from_${source}"

		# Make result directories 
		while IFS= read -r sub; do
			mkdir -p "$MINIMAP_DIR/$pair/$sub"
		done <<< "$TEMPLATE"

		mrna_fa=$(find  "$INPUT_DIR/$source" -name "*.mRNA.longest.exons.fa")
		lnc_fa=$(find   "$INPUT_DIR/$source" -name "*.lnc_RNA.longest.nonoverlapping.exons.fa")
		decoy_fa=$(find "$INPUT_DIR/$source" -name "*.decoy.longest.fa")

		ref=$(realpath "$target_fna")
		bam_pcg=$(realpath  -m "$MINIMAP_DIR/$pair/pcg/aligned.bam")
		bam_lnc=$(realpath  -m "$MINIMAP_DIR/$pair/lnc/aligned.bam")
		bam_decoy=$(realpath -m "$MINIMAP_DIR/$pair/decoy/aligned.bam")

		echo bash ./00_run_minimap_base.sh "$ref" "$(realpath "$mrna_fa")"  "$bam_pcg"
		echo bash ./00_run_minimap_base.sh "$ref" "$(realpath "$lnc_fa")"   "$bam_lnc"
		echo bash ./00_run_minimap_base.sh "$ref" "$(realpath "$decoy_fa")" "$bam_decoy"
	done
done
