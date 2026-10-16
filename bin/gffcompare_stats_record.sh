#!/usr/bin/env bash
# Write the `accuracy` data record (fomo-data/1) of a gffcompare .stats file: Sn/Pr per level,
# the query and reference counts, the matching counts and the missed/novel counts. Runs in the
# gffcompare image, which has no Python; fomo_stats.parse_gffcompare_stats is the same parser.
#
# Usage: gffcompare_stats_record.sh '<identity JSON>' <x.gffcompare.stats> > <record.data.json>

set -euo pipefail
export LC_ALL=C

if [[ $# -ne 2 ]]; then
    echo "Usage: $(basename "$0") '<identity JSON>' <x.gffcompare.stats> > <record.data.json>" >&2
    exit 1
fi

awk -v identity="$1" '
    function num(x) { return x ~ /^-?[0-9]+(\.[0-9]+)?$/ }
    function add(key, value) {
        sep = (counts == "") ? "" : ","
        counts = counts sep "\"" key "\":" value
    }
    /level:/ {
        line = $0; gsub(/\|/, "", line); sub(/Intron chain/, "Intron_chain", line)
        n = split(line, f, " ")
        if (n >= 4 && num(f[3]) && num(f[4])) {
            sep = (levels == "") ? "" : ","
            levels = levels sep "\"" f[1] "\":{\"sn\":" f[3] ",\"pr\":" f[4] "}"
        }
        next
    }
    /(Query|Reference) mRNAs :/ {
        line = $0; sub(/^#/, "", line); n = split(line, f, " ")
        key = (f[1] == "Query") ? "query" : "reference"
        add(key "_mrnas", f[4]); add(key "_loci", f[6])
        next
    }
    /^[ \t]*Matching (intron chains|transcripts|loci):/ {
        line = $0; sub(/^[ \t]*Matching /, "", line); split(line, f, ":")
        key = f[1]; gsub(/ /, "_", key)
        add("matching_" key, $NF)
        next
    }
    /^[ \t]*(Missed|Novel) (exons|introns|loci):/ {
        split($3, part, "/")
        key = tolower($1) "_" substr($2, 1, length($2) - 1)
        add(key, part[1]); add(key "_total", part[2])
    }
    END {
        printf "{\"schema\":\"fomo-data/1\",\"section\":\"accuracy\",\"identity\":%s,", identity
        printf "\"data\":{\"levels\":{%s},\"counts\":{%s}}}\n", levels, counts
    }
' "$2"
