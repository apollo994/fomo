#!/usr/bin/env bash
# nextflow lint over the whole pipeline, failing on warnings as well as errors.
set -euo pipefail
export LC_ALL=C

cd "$(dirname "$0")/.."
if ! out=$(nextflow lint -o concise main.nf workflows subworkflows modules conf nextflow.config 2>&1); then
    echo "$out"
    exit 1
fi
echo "$out"
if grep -q '^Warn ' <<< "$out"; then
    echo "nextflow lint: warnings are errors here" >&2
    exit 1
fi
