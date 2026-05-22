#!/bin/bash
#
# Emits one `01_convert_bam_to_gff.sh` invocation per BAM produced under
# ../../minimap_models/. The resulting GFF3 is written next to the BAM as
# <bam>.gff3 (so e.g. aligned.bam -> aligned.bam.gff3).
#
# Usage:
#   bash ./02_convert_bam_to_gff_generator.sh                       # print commands to stdout
#   bash ./02_convert_bam_to_gff_generator.sh > bam2gff_tasks.txt   # capture for slurm array
#   bash ./submit.sh bam2gff_tasks.txt bam2gff                      # launch via send_array.sh

set -euo pipefail

MINIMAP_DIR="../../minimap_models"

find "$MINIMAP_DIR" -type f -name "*.bam" | sort | while read -r bam; do
	bam_abs=$(realpath "$bam")
	gff_abs=$(realpath -m "${bam}.gff3")
	echo "bash ./01_convert_bam_to_gff.sh $bam_abs > $gff_abs"
done
