#!/usr/bin/env bash
# Compare every projected model set of one (target, track) with the target's reference: one
# gffcompare run per manifest row, and with --identity one `accuracy` data record per run.
#
# The manifest (gff_by_source.py --manifest) is a TSV `file<TAB>source`; the source is a species
# or the pooled allModels_raw. Each run reads the query as
# <target>.from_<source>.<track>.projected.gff3 (a symlink when the file has another name),
# because gffcompare names its .tmap/.refmap after the query and those files are published.
# Separate runs, not gffcompare's multi-query mode, which pools the queries into one Sn/Pr.
# Rationale: plans/18.
#
# Usage: gffcompare_batch.sh --manifest M --reference R --target T --track K
#                            [--identity '<JSON>'] [-- <gffcompare options>]

set -euo pipefail
export LC_ALL=C
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage="Usage: $(basename "$0") --manifest M --reference R --target T --track K [--identity JSON] [-- options]"
manifest="" reference="" target="" track="" identity=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)  manifest=$2; shift 2 ;;
        --reference) reference=$2; shift 2 ;;
        --target)    target=$2; shift 2 ;;
        --track)     track=$2; shift 2 ;;
        --identity)  identity=$2; shift 2 ;;
        --)          shift; break ;;
        *)           echo "$usage" >&2; exit 1 ;;
    esac
done
if [[ -z $manifest || -z $reference || -z $target || -z $track ]]; then
    echo "$usage" >&2
    exit 1
fi
if [[ $(head -n 1 "$manifest") != $'file\tsource' ]]; then
    echo "Error: $manifest is not a file/source manifest" >&2
    exit 1
fi

while IFS=$'\t' read -r file source; do
    prefix="${target}.from_${source}.${track}"
    query="${prefix}.projected.gff3"
    if [[ $file != "$query" && ! -e $query ]]; then
        ln -s "$file" "$query"
    fi
    gffcompare "$@" -r "$reference" -o "${prefix}.gffcompare" "$query"
    mv "${prefix}.gffcompare" "${prefix}.gffcompare.stats"
    if [[ -n $identity ]]; then
        "$here/gffcompare_stats_record.sh" "{\"source\":\"${source}\",${identity:1}" "${prefix}.gffcompare.stats" \
            > "${prefix}.accuracy.data.json"
    fi
done < <(tail -n +2 "$manifest")
