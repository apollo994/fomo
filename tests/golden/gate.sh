#!/usr/bin/env bash
# G4 on the CRG cluster: run V1-V4 as SLURM jobs (4 cpu each) and compare with the golden snapshots.
# Usage: tests/golden/gate.sh SCRATCH_DIR   (blocks until all four finish; exit 1 on any difference)
set -euo pipefail
export LC_ALL=C
scratch=$(realpath -m "$1")
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$scratch/logs"
jobs=()
for v in V1 V2 V3 V4; do
    rm -rf "${scratch:?}/$v"
    jobs+=("$(sbatch --parsable --qos=normal --partition=genoa64 --cpus-per-task=4 --mem=16G --time=2:00:00 \
        -J "gate_$v" -o "$scratch/logs/$v.out" -e "$scratch/logs/$v.err" \
        --wrap "$here/run_variant.sh $v $scratch/$v" 2> /dev/null)")
done
while squeue -h -j "$(IFS=,; echo "${jobs[*]}")" 2> /dev/null | grep -q .; do sleep 15; done
status=0
for v in V1 V2 V3 V4; do
    if [[ ! -s $scratch/$v/snapshot.tsv ]]; then
        echo "$v: no snapshot (run failed, see $scratch/$v/nextflow.stdout)"; status=1; continue
    fi
    if python3 "$here/compare.py" "$here/$v.tsv" "$scratch/$v/snapshot.tsv" > "$scratch/$v/compare.txt" 2>&1; then
        echo "$v: identical"
    else
        echo "$v: DIFFERENT"; cat "$scratch/$v/compare.txt"; status=1
    fi
done
exit $status
