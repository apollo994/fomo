#!/usr/bin/env bash
# Write the `alignment` data record (fomo-data/1) of a `samtools stats` file: the SN fields of
# the alignment table and reads mapped %. Runs in the samtools image, which has no Python.
#
# Usage: samtools_stats_record.sh '<identity JSON>' <file.stats> > <record.data.json>

set -euo pipefail
export LC_ALL=C

if [[ $# -ne 2 ]]; then
    echo "Usage: $(basename "$0") '<identity JSON>' <file.stats> > <record.data.json>" >&2
    exit 1
fi

awk -F'\t' -v identity="$1" '
    function value(key) { return (key in sn) ? sn[key] : "null" }
    $1 == "SN" {
        key = $2; sub(/:$/, "", key)
        val = $3; sub(/[ \t]*#.*$/, "", val); gsub(/^[ \t]+|[ \t]+$/, "", val)
        sn[key] = val
    }
    END {
        pct = "null"
        if (("raw total sequences" in sn) && sn["raw total sequences"] + 0 > 0)
            pct = sprintf("%.2f", 100 * sn["reads mapped"] / sn["raw total sequences"])
        printf "{\"schema\":\"fomo-data/1\",\"section\":\"alignment\",\"identity\":%s,\"data\":{", identity
        printf "\"raw_total_sequences\":%s,\"reads_mapped\":%s,\"reads_mapped_percent\":%s,", value("raw total sequences"), value("reads mapped"), pct
        printf "\"non_primary_alignments\":%s,\"error_rate\":%s,", value("non-primary alignments"), value("error rate")
        printf "\"average_quality\":%s,\"average_length\":%s}}\n", value("average quality"), value("average length")
    }
' "$2"
