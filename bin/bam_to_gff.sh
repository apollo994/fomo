#!/usr/bin/env bash
# Convert a minimap2 spliced-alignment BAM into a GFF3 of projected transcripts.
#
# Read names are expected in the form  <tid>|<gene_class>|<source>
# (produced by RENAME_FASTA_HEADERS). The transcript ID in the GFF is
# <tid>|<source> — NOT the bare <tid> alone: source transcript ids are not
# guaranteed unique across species (seen for real: two Drosophila species both
# used "lnc_RNA1412" as their own id), and every source's models land in one
# merged GFF3 downstream (plans/15), so bare tid would collide. <source> and
# <gene_class> are also promoted to 9th-column attributes as before. NM, AS
# and de tags from the BAM are also copied onto the transcript row.
#
# By default, single-exon projected models (alignments with no N/splice gap,
# i.e. one exon block) are discarded — the source annotations are spliced
# transcripts, so a single-exon projection usually reflects a collapsed or
# spurious mapping. Pass --include-single-exon to keep them.
#
# Usage: bam_to_gff.sh [--include-single-exon] <input.bam> > <output.gff3>

set -euo pipefail

include_single_exon=0
bam_file=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --include-single-exon)
            include_single_exon=1
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [--include-single-exon] <input.bam> > <output.gff3>" >&2
            exit 0
            ;;
        -*)
            echo "Error: unknown option '$1'" >&2
            echo "Usage: $0 [--include-single-exon] <input.bam>" >&2
            exit 1
            ;;
        *)
            if [[ -n "$bam_file" ]]; then
                echo "Error: unexpected extra argument '$1'" >&2
                exit 1
            fi
            bam_file="$1"
            shift
            ;;
    esac
done

if [[ -z "$bam_file" ]]; then
    echo "Usage: $0 [--include-single-exon] <input.bam>" >&2
    exit 1
fi
if [[ ! -f "$bam_file" ]]; then
    echo "Error: input BAM '$bam_file' not found." >&2
    exit 1
fi

printf '##gff-version 3\n'

# Single pass over the BAM: drop unmapped (0x4), secondary (0x100) and
# supplementary (0x800) records (-F 2308), then walk the CIGAR to derive
# exon blocks in target genomic coordinates.
samtools view -F 2308 "$bam_file" |
awk -v OFS='\t' -v include_single_exon="$include_single_exon" '
    {
        qname  = $1
        flag   = $2
        rname  = $3
        pos    = $4 + 0          # 1-based leftmost mapping pos
        mapq   = $5
        cigar  = $6
        strand = (int(flag / 16) % 2 == 1 ? "-" : "+")

        # ---- read-name → (tid, gene_class, source) -------------------------
        np = split(qname, parts, "|")
        tid    = parts[1]
        gclass = (np >= 2 ? parts[2] : ".")
        src    = (np >= 3 ? parts[3] : ".")
        if (tid ~ /[;=]/) {
            print "Error: tid contains ; or = which would corrupt GFF attributes: " tid > "/dev/stderr"
            exit 1
        }

        # ---- pull NM, AS, de from optional SAM tags ------------------------
        nm = "."; as = "."; de = "."
        for (i = 12; i <= NF; i++) {
            if      ($i ~ /^NM:i:/) nm = substr($i, 6)
            else if ($i ~ /^AS:i:/) as = substr($i, 6)
            else if ($i ~ /^de:f:/) de = substr($i, 6)
        }

        # ---- walk CIGAR to derive exon blocks ------------------------------
        # Block accumulator in genomic (reference) coordinates, plus three
        # running totals used to derive the size/coverage attributes:
        #   query_len     = ∑ M, I, S, =, X       (read length)
        #   aligned_bases = ∑ M, =, X             (query bases matched on ref)
        #   projected_len = ∑ M, D, =, X          (exon-block bases on ref)
        block_start   = 0
        block_len     = 0
        n_blocks      = 0
        ref_pos       = pos
        query_len     = 0
        aligned_bases = 0
        projected_len = 0
        rest          = cigar
        while (match(rest, /^[0-9]+[MIDNSHP=X]/) > 0) {
            tok = substr(rest, 1, RLENGTH)
            op  = substr(tok, RLENGTH, 1)
            len = substr(tok, 1, RLENGTH - 1) + 0
            rest = substr(rest, RLENGTH + 1)

            if (op == "M" || op == "=" || op == "X") {
                if (block_len == 0) block_start = ref_pos
                block_len     += len
                ref_pos       += len
                query_len     += len
                aligned_bases += len
                projected_len += len
            }
            else if (op == "D") {
                if (block_len == 0) block_start = ref_pos
                block_len     += len
                ref_pos       += len
                projected_len += len   # consumes reference only
            }
            else if (op == "N") {
                if (block_len > 0) {
                    block_ends[++n_blocks] = block_start "," (block_start + block_len - 1)
                }
                block_start = 0
                block_len   = 0
                ref_pos    += len
            }
            else if (op == "I" || op == "S") {
                query_len += len       # consumes query only
            }
            else if (op == "H" || op == "P") {
                # consumes neither
            }
        }
        if (block_len > 0) {
            block_ends[++n_blocks] = block_start "," (block_start + block_len - 1)
        }
        if (n_blocks == 0) next   # nothing to emit

        # Drop single-exon models unless explicitly kept. Source transcripts are
        # spliced, so a single exon block (no N gap) is usually a collapsed or
        # spurious projection.
        if (n_blocks == 1 && include_single_exon == 0) next

        # transcript/gene span = first block start .. last block end
        split(block_ends[1],        first, ",")
        split(block_ends[n_blocks], last,  ",")
        tx_start = first[1] + 0
        tx_end   = last[2]  + 0

        pct_aligned   = (query_len > 0 ? sprintf("%.1f", 100.0 * aligned_bases / query_len) : ".")
        pct_projected = (query_len > 0 ? sprintf("%.1f", 100.0 * projected_len / query_len) : ".")

        # Globally unique across every source merged into allModels — see the
        # header comment. The duplicate-ID guard in gff_by_source.py is what
        # would catch a regression here.
        uid = tid "|" src

        gene_attrs = "ID=gene-" uid ";source=" src ";gene_class=" gclass
        tx_attrs   = "ID=" uid ";Parent=gene-" uid \
                     ";source=" src ";gene_class=" gclass \
                     ";query_length=" query_len \
                     ";aligned_bases=" aligned_bases \
                     ";pct_aligned=" pct_aligned \
                     ";projected_length=" projected_len \
                     ";pct_projected=" pct_projected \
                     ";NM=" nm ";AS=" as ";de=" de

        print rname, "fomo", "gene",       tx_start, tx_end, mapq, strand, ".", gene_attrs
        print rname, "fomo", "transcript", tx_start, tx_end, mapq, strand, ".", tx_attrs

        for (i = 1; i <= n_blocks; i++) {
            split(block_ends[i], b, ",")
            ex_attrs = "ID=" uid ".exon" i ";Parent=" uid
            print rname, "fomo", "exon", b[1] + 0, b[2] + 0, mapq, strand, ".", ex_attrs
        }

        delete block_ends
    }
'
