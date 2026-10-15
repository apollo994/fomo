#!/usr/bin/env bash
# Reference filter on BOTH strands (curate_models.py --strand both) vs the default same-strand
# filter, at the default --min-species 2, for every benchmarked target's lncRNA track.
# Usage: strand_both.sh RESULTS_DIR OUT_DIR SAMPLESHEET
set -euo pipefail
RES=$1; OUT=$2; SHEET=$3
REPO=/users/rg/fzanarello/pipelines/fomo
IMG=/software/rg/images/community.wave.seqera.io-library-python_gffcompare_gffread-970e511c83295d93.img
mkdir -p "$OUT"
tail -n +2 "$SHEET" | while IFS=, read -r sp role _fa gff; do
    [[ $role == both ]] || continue
    echo "$sp,$gff"
done | xargs -P "${SLURM_CPUS_PER_TASK:-4}" -I{} bash -c '
    sp=${1%%,*}; gff=${1#*,}
    singularity exec -B /nfs,/users '"$IMG"' python3 '"$REPO"'/bin/curate_models.py \
        --gff '"$RES"'/targets/$sp/annotation/$sp.allModels.lnc_RNA.raw.gff3 --ref "$gff" \
        --prefix '"$OUT"'/$sp.lnc_RNA.strand_both --target $sp --track lnc_RNA \
        --min-species 2 --strand both --skip-gffread > '"$OUT"'/$sp.log 2>&1
    rm -f '"$OUT"'/$sp.lnc_RNA.strand_both.gff3.gz '"$OUT"'/$sp.lnc_RNA.strand_both.merged.gff3.gz' _ {}
