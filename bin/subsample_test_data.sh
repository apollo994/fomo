#!/usr/bin/env bash
# One-off script: subsample test data to the longest contig per species.
# Replaces symlinks in assets/test_data/ with committed subsampled files.
# Run once from the repo root: bash bin/subsample_test_data.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ANNO_DIR="$REPO_DIR/annotation_downloads"
TEST_DATA_DIR="$REPO_DIR/assets/test_data"

# "species:assembly" pairs
PAIRS=(
    "Lycaena_phlaeas_282391:GCA_905333005.2"
    "Lycaena_alciphron_282377:GCA_964273875.1"
    "Lycaena_hippothoe_580924:GCA_964264015.1"
    "Lycaena_thersamon_265360:GCA_965194355.1"
    "Lycaena_virgaureae_282395:GCA_964263885.1"
)

for pair in "${PAIRS[@]}"; do
    species="${pair%%:*}"
    acc="${pair##*:}"
    src_dir="$ANNO_DIR/$species/$acc"
    dst_dir="$TEST_DATA_DIR/$species/$acc"

    echo "=== $species ($acc) ==="

    fna_gz=$(find "$src_dir" -name "*.fna.gz" | head -1)
    gff3_gz=$(find "$src_dir" -name "*.aliasMatch.gff3.gz" | head -1)

    [[ -z "$fna_gz"  ]] && { echo "ERROR: no .fna.gz in $src_dir" >&2; exit 1; }
    [[ -z "$gff3_gz" ]] && { echo "ERROR: no .gff3.gz in $src_dir" >&2; exit 1; }

    # Replace symlink with real directory
    if [[ -L "$TEST_DATA_DIR/$species" ]]; then
        rm "$TEST_DATA_DIR/$species"
    fi
    mkdir -p "$dst_dir"

    out_fna="$dst_dir/${species}.fna.gz"
    out_gff3="$dst_dir/${species}.gff3.gz"

    # Find the longest contig name (two-pass-free: accumulate lengths then pick max)
    echo "  Finding longest contig ..."
    LONGEST=$(gunzip -c "$fna_gz" | awk '
        /^>/ {
            if (name != "") lens[name] = len
            name = substr($1, 2)
            len = 0
            next
        }
        { len += length($0) }
        END {
            if (name != "") lens[name] = len
            max_len = 0; best = ""
            for (n in lens) {
                if (lens[n] > max_len) { max_len = lens[n]; best = n }
            }
            print best
        }
    ')
    echo "  Longest contig: $LONGEST"

    # Extract the single contig FASTA
    echo "  Writing $out_fna ..."
    gunzip -c "$fna_gz" | awk -v seq="$LONGEST" '
        /^>/ { found = (substr($1, 2) == seq) }
        found { print }
    ' | bgzip > "$out_fna"

    # Subsample GFF3 to matching contig
    echo "  Writing $out_gff3 ..."
    gunzip -c "$gff3_gz" \
        | awk -v c="$LONGEST" -F'\t' 'BEGIN{OFS="\t"} /^#/ || $1==c' \
        | gzip > "$out_gff3"

    echo "  Done."
done

echo ""
echo "Subsampling complete."
echo "Next: commit assets/test_data/, then run: nextflow run . -profile test,docker"
