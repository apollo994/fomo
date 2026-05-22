#!/bin/bash
# Builds the annotations.csv samplesheet for the TPFP validateSJ run.
# References every <basename>.{TP,FP}.gff3 in benchmark/TPFP_analysis/input/
# with annotation_id = <basename>.{TP|FP}.
set -euo pipefail

ROOT="/users/project/phyloENCODE_009142_no_backup/general_analyses/crossmapping/lycaenidae_benchmark"
INPUT_DIR="$ROOT/benchmark/TPFP_analysis/input"
OUT="$ROOT/benchmark/TPFP_analysis/validateJS/annotations.csv"

mkdir -p "$(dirname "$OUT")"

{
	echo "annotation_id,gff3"
	for kind in TP FP
	do
		for gff in "$INPUT_DIR"/*."$kind".gff3
		do
			[[ -e "$gff" ]] || continue
			base=$(basename "$gff" .gff3)        # e.g. liftoff_..._pcg.TP
			echo "$base,$(realpath "$gff")"
		done
	done
} > "$OUT"

echo "wrote $OUT ($(($(wc -l < "$OUT") - 1)) annotations)"
