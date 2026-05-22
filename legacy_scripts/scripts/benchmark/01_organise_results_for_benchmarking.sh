#!/bin/bash
# Collect minimap/liftoff GFF3 outputs into benchmark/input/ as symlinks named
#   {tool}_{target}_from_{source}_{genetype}.gff3
#
# Source layouts:
#   minimap_models/<target>_from_<source>/<genetype>/aligned.bam.gff3
#   liftoff_models/<target>_from_<source>/<genetype>/lifted.<genetype>.gff3
# Genetypes: pcg, lnc, decoy
set -euo pipefail

ROOT="/users/project/phyloENCODE_009142_no_backup/general_analyses/crossmapping/lycaenidae_benchmark"
INPUT_DIR="$ROOT/benchmark/input"
GENETYPES=(pcg lnc decoy)

mkdir -p "$INPUT_DIR"

n_linked=0
n_missing=0

for tool in minimap liftoff; do
    src_root="$ROOT/${tool}_models"
    [[ -d "$src_root" ]] || { echo "skip: $src_root not found"; continue; }

    for model_dir in "$src_root"/*/; do
        model="$(basename "$model_dir")"
        # only process directories matching <target>_from_<source>
        [[ "$model" == *_from_* ]] || continue

        for gt in "${GENETYPES[@]}"; do
            case "$tool" in
                minimap) src="$model_dir$gt/aligned.bam.gff3" ;;
                liftoff) src="$model_dir$gt/lifted.${gt}.gff3" ;;
            esac

            link="$INPUT_DIR/${tool}_${model}_${gt}.gff3"

            if [[ -s "$src" ]]; then
                ln -sfn "$src" "$link"
                n_linked=$((n_linked + 1))
            else
                echo "missing: $src"
                n_missing=$((n_missing + 1))
            fi
        done
    done
done

echo
echo "Linked  : $n_linked"
echo "Missing : $n_missing"
echo "Output  : $INPUT_DIR"
