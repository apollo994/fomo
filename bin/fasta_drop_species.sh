#!/usr/bin/env bash
# Write a FASTA without the records of one species: those whose header id (`>tid|type|species`,
# RENAME_FASTA_HEADERS) has <species> as its last |-field.
#
# Usage: fasta_drop_species.sh <species> <in.fasta> > <out.fasta>

set -euo pipefail
export LC_ALL=C

if [[ $# -ne 2 ]]; then
    echo "Usage: $(basename "$0") <species> <in.fasta> > <out.fasta>" >&2
    exit 1
fi

awk -v s="$1" '/^>/ { n = split(substr($1, 2), f, "|"); keep = (f[n] != s) } keep' "$2"
