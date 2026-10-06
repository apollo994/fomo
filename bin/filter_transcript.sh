#!/usr/bin/env bash
# filter_transcript.sh
# Usage: filter_transcript.sh [--longest] [--counts FILE --sample NAME] \
#            input.gff[3] TRANSCRIPT_TYPE > output.gff
#
# Emits only **stranded, spliced (multi-exon) transcripts** of type
# TRANSCRIPT_TYPE along with their parent gene and child exon/CDS/UTR features.
#
# Stranded: the transcript's own strand is + or -, every child (exon/CDS/UTR)
# is on that same strand, and no parent has an undefined strand. Anything else
# is dropped. These are trans-spliced models — e.g. RefSeq plant mitochondrial
# nad1/nad2/nad5, whose exons sit on both strands and whose mRNA line carries
# strand "?". gffread -w hard-errors on "?" ("Error parsing strand"), and a
# mixed-strand model cannot be extracted as one contiguous cDNA nor projected
# as one alignment, so there is nothing useful to transfer. Seen for real on
# the Apiaceae run (Daucus carota, RefSeq GCF_001625215.2).
#
# Transcripts with fewer than 2 exons are dropped, which keeps source projections and the
# benchmarking reference apples-to-apples (gffcompare runs with -M
# downstream, and decoys already require >=2 exons in relocate_loci.py).
#
# --longest keeps ONE transcript per gene (plans/20_longest_isoform.md): among
# the gene's SPLICED isoforms of TRANSCRIPT_TYPE, the one with the greatest
# summed exon length — the spliced length, i.e. exactly what `gffread -w`
# extracts and minimap2 aligns (not the intron-inflated genomic span, and not
# CDS length, which would rank mRNA differently from lncRNA). Rules:
#   * selection runs AFTER the multi-exon filter, so a gene whose longest
#     isoform is mono-exonic keeps its longest *spliced* one instead of vanishing;
#   * a "gene" is the transcript's full Parent= value; a transcript with no
#     Parent competes only with itself;
#   * ties go to the lexicographically smallest transcript ID (byte order,
#     LC_ALL=C) — awk's for-in order is unspecified, so without an explicit rule
#     the pick could change between runs.
# Without --longest every spliced transcript is kept, exactly as before.
#
# --counts FILE --sample NAME append a one-row funnel table (MultiQC custom
# content, *_transcript_filter_mqc.tsv) with the counts at each step:
#   n_transcripts  candidates of TRANSCRIPT_TYPE
#   n_stranded     of those, consistently stranded (see above)
#   n_spliced      of those, with >= 2 exons
#   n_selected     emitted (== n_spliced without --longest)
set -euo pipefail
export LC_ALL=C

usage="Usage: $0 [--longest] [--counts FILE --sample NAME] input.gff[3] TRANSCRIPT_TYPE > output.gff"

longest=0
counts=""
sample=""
positional=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --longest)
      longest=1
      shift
      ;;
    --counts)
      counts="${2:-}"
      shift 2
      ;;
    --sample)
      sample="${2:-}"
      shift 2
      ;;
    -h|--help)
      echo "$usage" >&2
      exit 0
      ;;
    -*)
      echo "Error: unknown option '$1'" >&2
      echo "$usage" >&2
      exit 1
      ;;
    *)
      positional+=("$1")
      shift
      ;;
  esac
done

gff="${positional[0]:-}"
ttype="${positional[1]:-}"

if [[ -z "${gff}" || -z "${ttype}" || ${#positional[@]} -ne 2 ]]; then
  echo "$usage" >&2
  exit 1
fi
if [[ -n "${counts}" && -z "${sample}" ]]; then
  echo "Error: --counts requires --sample" >&2
  exit 1
fi

awk -F'\t' -v TTYPE="$ttype" -v LONGEST="$longest" -v COUNTS="$counts" -v SAMPLE="$sample" '
  BEGIN { OFS="\t" }

  # ---------- PASS 1: candidates, exon count/length and strands per parent ----------
  FNR==NR {
    if ($0 ~ /^#/ || NF < 9) next
    id=""; parent=""
    n=split($9,a,";")
    for (i=1;i<=n;i++) {
      if (a[i] ~ /^ID=/ && id == "")              { id=a[i]; sub(/^ID=/,"",id) }
      else if (a[i] ~ /^Parent=/ && parent == "") { parent=a[i]; sub(/^Parent=/,"",parent) }
    }
    # Any feature with an undefined strand ("?", ".") taints every transcript that
    # names it as a parent (gffread would choke on the parent line too).
    if (id != "" && $7 != "+" && $7 != "-") bad_strand_id[id] = 1
    if ($3 == TTYPE) {
      if (id != "") { tx_candidate[id] = parent; tx_strand[id] = $7 }
    }
    else if (parent != "") {
      # Record every child strand per parent (trans-spliced = children on both
      # strands), and tally exons (count + summed length) for the multi-exon filter
      # and the longest-isoform selection. An exon shared by several transcripts
      # (Parent=a,b) credits each of them.
      m=split(parent,p,",")
      for (j=1;j<=m;j++) {
        if (p[j]=="") continue
        if      ($7 == "+") child_plus[p[j]]  = 1
        else if ($7 == "-") child_minus[p[j]] = 1
        else                child_other[p[j]] = 1
        if ($3 == "exon") { exon_count[p[j]]++; exon_len[p[j]] += $5 - $4 + 1 }
      }
    }
    next
  }

  # ---------- Boundary: first line of second file copy ----------
  # Promote multi-exon candidates (with --longest: only the longest one per gene)
  # into the final tx[] / parent_of_tx[] sets.
  FNR==1 && !promoted {
    promoted=1
    n_transcripts=0; n_stranded=0; n_spliced=0; n_selected=0
    for (cid in tx_candidate) {
      n_transcripts++
      st = tx_strand[cid]
      if (st != "+" && st != "-") continue
      if ((cid in child_other) || (st == "+" && (cid in child_minus)) \
                               || (st == "-" && (cid in child_plus))) continue
      bad_parent = 0
      m=split(tx_candidate[cid],p,",")
      for (j=1;j<=m;j++) if (p[j] != "" && (p[j] in bad_strand_id)) bad_parent = 1
      if (bad_parent) continue
      n_stranded++
      if (exon_count[cid] < 2) continue
      n_spliced++
      if (!LONGEST) { tx[cid] = 1; continue }
      key = (tx_candidate[cid] != "") ? tx_candidate[cid] : cid
      if (!(key in best) \
          || exon_len[cid] > best_len[key] \
          || (exon_len[cid] == best_len[key] && cid < best[key])) {
        best[key] = cid
        best_len[key] = exon_len[cid]
      }
    }
    if (LONGEST) for (key in best) tx[best[key]] = 1
    for (cid in tx) {
      n_selected++
      cparent = tx_candidate[cid]
      if (cparent != "") {
        m=split(cparent,p,",")
        for (j=1;j<=m;j++) if (p[j]!="") parent_of_tx[p[j]] = 1
      }
    }
    if (COUNTS != "") {
      print "Sample", "feature", "n_transcripts", "n_stranded", "n_spliced", "n_selected" > COUNTS
      print SAMPLE, TTYPE, n_transcripts, n_stranded, n_spliced, n_selected > COUNTS
      close(COUNTS)
    }
  }

  # ---------- helpers ----------
  function get_attr(field, key,   n,i,a,kv) {
    n=split(field,a,";")
    for(i=1;i<=n;i++){
      if (a[i] ~ ("^" key "=")) {
        kv=a[i]
        sub("^" key "=", "", kv)
        return kv
      }
    }
    return ""
  }

  # ---------- PASS 2: print lines matching your rules ----------
  /^#/ { next }
  NF < 9 { next }

  {
    id     = get_attr($9, "ID")
    parent = get_attr($9, "Parent")

    # Rule 1: print if this feature is a parent of a selected transcript
    # (works for gene, nc_gene, or any other parent feature type)
    if (id != "" && parent_of_tx[id]) { print; next }

    # Rule 2: print if this feature ID is one of the selected transcripts (the ttype lines)
    if (id != "" && tx[id]) { print; next }

    # Rule 3: print if this feature is a child of a selected transcript (exon, CDS, UTR, etc.)
    if (parent != "") {
      m=split(parent,p,",")
      for (j=1;j<=m;j++) {
        if (tx[p[j]]) { print; next }
      }
    }
  }

  # An empty input never reaches the boundary block; still emit the funnel row so
  # the module output exists (and honestly reports zero).
  END {
    if (!promoted && COUNTS != "") {
      print "Sample", "feature", "n_transcripts", "n_stranded", "n_spliced", "n_selected" > COUNTS
      print SAMPLE, TTYPE, 0, 0, 0, 0 > COUNTS
    }
  }
' "$gff" "$gff"
