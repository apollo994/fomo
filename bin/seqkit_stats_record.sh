#!/usr/bin/env bash
# Write the `fasta_stats` data record (fomo-data/1) of a `seqkit stats --tabular` file: every
# column but `file`, numbers as numbers. Runs in the seqkit image, which has no Python.
#
# Usage: seqkit_stats_record.sh '<identity JSON>' <stats.tsv> > <record.data.json>

set -euo pipefail
export LC_ALL=C

if [[ $# -ne 2 ]]; then
    echo "Usage: $(basename "$0") '<identity JSON>' <stats.tsv> > <record.data.json>" >&2
    exit 1
fi

awk -F'\t' -v identity="$1" '
    NR == 1 { for (i = 1; i <= NF; i++) col[i] = $i; next }
    NR == 2 {
        out = ""
        for (i = 1; i <= NF; i++) {
            if (col[i] == "file") continue
            v = ($i ~ /^-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?$/) ? $i : "\"" $i "\""
            sep = (out == "") ? "" : ","
            out = out sep "\"" col[i] "\":" v
        }
    }
    END {
        printf "{\"schema\":\"fomo-data/1\",\"section\":\"fasta_stats\",\"identity\":%s,\"data\":{%s}}\n", identity, out
    }
' "$2"
