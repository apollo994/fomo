#!/bin/bash
# Build ENA fastq.gz URLs from accession IDs.
# Reads:  RNAseq/Polyommatus_icarus_265386_RNAseq.txt
# Writes: RNAseq/urls.txt   (one URL per line, R1 then R2 for each accession)
#
# Pattern for 7-digit ERR/SRR accessions:
#   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/<first6>/00<last>/<acc>/<acc>_{1,2}.fastq.gz
set -euo pipefail

ROOT="/users/project/phyloENCODE_009142_no_backup/general_analyses/crossmapping/lycaenidae_benchmark"
ACC_LIST="$ROOT/RNAseq/Polyommatus_icarus_265386_RNAseq.txt"
OUT_DIR="$ROOT/RNAseq"

while read -r acc; do
    acc="${acc//[$'\r\n\t ']}"          # strip whitespace / CR
    [[ -z "$acc" ]] && continue
    prefix="${acc:0:6}"                  # e.g. ERR316
    last="${acc: -1}"                    # last digit
    base="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/00${last}/${acc}"
    echo "wget -c -P ${OUT_DIR} ${base}/${acc}_1.fastq.gz"
    echo "wget -c -P ${OUT_DIR} ${base}/${acc}_2.fastq.gz"
done < "$ACC_LIST"

