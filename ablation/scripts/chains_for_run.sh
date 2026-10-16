#!/usr/bin/env bash
# Re-curate every benchmarked target of a FOMO run at --min-species 1 and dump the chain audit
# table (--chains-tsv), with and without the reference filter. One run per track gives every
# supported chain with its supporting species, so any support threshold and any donor subset
# can be evaluated afterwards without re-running curation (ablation/scripts/ablate.py).
#
# Runs bin/curate_models.py inside the pipeline's own CURATE_MODELS container, so gffcompare
# (0.12.6) and python (3.11) are the versions the pipeline uses.
#
# Usage: chains_for_run.sh RESULTS_DIR OUT_DIR SAMPLESHEET "TRACKS" [JOBS]
#   TRACKS e.g. "lnc_RNA lnc_RNA.decoy mRNA"
set -euo pipefail

RESULTS=$1; OUT=$2; SHEET=$3; TRACKS=$4; JOBS=${5:-${SLURM_CPUS_PER_TASK:-4}}
REPO=${REPO:-/users/rg/fzanarello/pipelines/fomo}
IMG=${IMG:-/software/rg/images/community.wave.seqera.io-library-python_gffcompare_gffread-970e511c83295d93.img}
mkdir -p "$OUT"

# One line per job: target \t track \t mode \t ref
jobs=$OUT/jobs.tsv
: > "$jobs"
tail -n +2 "$SHEET" | while IFS=, read -r sp role _fa gff; do
    [[ $role == both ]] || continue
    for tr in $TRACKS; do
        [[ -s $RESULTS/targets/$sp/annotation/$sp.allModels.$tr.raw.gff3 ]] || { echo "missing $sp $tr" >&2; continue; }
        printf '%s\t%s\tnoref\t-\n' "$sp" "$tr" >> "$jobs"
        printf '%s\t%s\tref\t%s\n'  "$sp" "$tr" "$gff" >> "$jobs"
    done
done

run_one() {
    IFS=$'\t' read -r sp tr mode ref <<< "$1"
    local pre=$OUT/$sp/$sp.$tr.$mode
    [[ -s $pre.chains.tsv.gz ]] && return 0
    mkdir -p "$OUT/$sp"
    local refarg=()
    [[ $mode == ref ]] && refarg=(--ref "$ref")
    singularity exec -B /nfs,/users "$IMG" python3 "$REPO/bin/curate_models.py" \
        --gff "$RESULTS/targets/$sp/annotation/$sp.allModels.$tr.raw.gff3" \
        --prefix "$pre" --target "$sp" --track "$tr" \
        --min-species 1 --chains-tsv --skip-gffread "${refarg[@]}" \
        > "$pre.log" 2>&1 || { echo "FAILED $sp $tr $mode" >&2; return 1; }
    # Only the audit table and the report are needed; the curated GFF3s are bulky at min 1.
    rm -f "$pre.gff3.gz" "$pre.merged.gff3.gz"
}
export -f run_one
export OUT RESULTS REPO IMG

tr '\n' '\0' < "$jobs" | xargs -0 -P "$JOBS" -I{} bash -c 'run_one "$1"' _ {}
echo "done: $(wc -l < "$jobs") jobs"
