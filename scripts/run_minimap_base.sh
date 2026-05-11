#!/bin/bash

# Exit immediately on error, undefined variable, or pipe failure
set -euo pipefail

# Load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate isoquant

# Check inputs
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <reference.fa> <reads.fastq> <output.bam>"
    exit 1
fi

REF="$1"
READS="$2"
BAM="$3"
SAM="${BAM%.bam}.sam"

# Validate input files
if [[ ! -f "$REF" ]]; then
    echo "Error: Reference genome '$REF' not found."
    exit 1
fi

if [[ ! -f "$READS" ]]; then
    echo "Error: Reads file '$READS' not found."
    exit 1
fi

# Check if output already exists
if [[ -f "$SAM" || -f "$BAM" || -f "${BAM}.bai" ]]; then
    echo "Output files already exist. Please remove or rename them first."
    exit 1
fi

# Log start time
echo "Starting minimap2 alignment at $(date)"
echo "Reference: $REF"
echo "Reads: $READS"
echo "Output BAM: $BAM"

# Run minimap2
time minimap2 -ax splice -t 4 "$REF" "$READS" > "$SAM"

# Sort and compress
echo "Sorting and compressing SAM to BAM at $(date)"
time samtools sort -@ 1 -o "$BAM" "$SAM"

# Index BAM
echo "Indexing BAM at $(date)"
time samtools index -@ 1 "$BAM"

# Get stats for all seq
echo "Extracting stats $(date)"
time samtools stats -@ 1 "$BAM" > "$BAM".stats.all
grep ^SN "$BAM".stats.all | cut -f 2- > "$BAM".stats.all.summary

# Get stats for primary mappings only
echo "Extracting stats $(date)"
time samtools stats -F 2308 -@ 1 "$BAM" > "$BAM".stats.primary
grep ^SN "$BAM".stats.primary | cut -f 2- > "$BAM".stats.primary.summary

# Remove intermediate SAM file
echo "Removing intermediate SAM file"
rm "$SAM"

# Log end time
echo "All steps completed at $(date)"
