#!/usr/bin/env bash
# Run one golden variant of the test profile and snapshot its deliverables (plans/27_refactor.md §5).
#
# Usage: run_variant.sh VARIANT OUTDIR [PROFILE]
#   VARIANT  V1 default | V2 --include_mrna --include_decoy | V3 reference-less target |
#            V4 = V1 with plain-gzip FASTAs and uncompressed GFF3s (must equal V1)
#   PROFILE  default test,singularity (CI: test,docker)
# Writes OUTDIR/results, OUTDIR/snapshot.tsv and OUTDIR/snapshot.info.tsv. Inside a SLURM job the
# resource limits follow the allocation.
set -euo pipefail
export LC_ALL=C

variant=$1
outdir=$(realpath -m "$2")
profile=${3:-test,singularity}
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)

sheet=assets/samplesheet.csv
extra=()
case "$variant" in
    V1 | V4) ;;
    V2) extra=(--include_mrna --include_decoy) ;;
    V3) sheet=assets/samplesheet_target_noref.csv ;;
    *) echo "unknown variant: $variant" >&2; exit 1 ;;
esac

mkdir -p "$outdir"
cd "$outdir"
sed "s#,assets/#,$repo/assets/#g" "$repo/$sheet" > samplesheet.csv

if [[ $variant == V4 ]]; then
    # plain gzip FASTA (not bgzip) and uncompressed GFF3: the production compression paths
    mkdir -p inputs
    {
        IFS= read -r header
        echo "$header"
        while IFS=, read -r species role fasta gff3; do
            fa=inputs/$(basename "$fasta")
            gzip -cd "$fasta" | gzip -n > "$fa"
            gff=""
            if [[ -n $gff3 ]]; then
                gff=inputs/$(basename "$gff3" .gz)
                gzip -cd "$gff3" > "$gff"
                gff=$outdir/$gff
            fi
            echo "$species,$role,$outdir/$fa,$gff"
        done
    } < samplesheet.csv > samplesheet.v4.csv
    mv samplesheet.v4.csv samplesheet.csv
fi

cfg=()
if [[ -n ${SLURM_CPUS_PER_TASK:-} ]]; then
    printf 'process.resourceLimits = [cpus: %s, memory: 16.GB, time: 2.h]\n' "$SLURM_CPUS_PER_TASK" > limits.config
    cfg=(-c limits.config)
fi

nextflow "${cfg[@]}" run "$repo" -profile "$profile" \
    --input samplesheet.csv --outdir results -w work "${extra[@]}" > nextflow.stdout 2>&1

python3 "$here/snapshot.py" results --info snapshot.info.tsv > snapshot.tsv
echo "snapshot: $outdir/snapshot.tsv ($(($(wc -l < snapshot.tsv) - 1)) deliverables)"
