#!/usr/bin/env python3
"""
Filter coding-potential transcripts out of a lncRNA FASTA using TD2 output.

A transcript is called CODING (and dropped) when it has at least one ORF in the
TD2.Predict `.pep` file that is BOTH:
  * complete (ORF type contains "complete"), and
  * PSAURON-confident (psauron_score >= --psauron-min).

TD2.Predict `.pep` headers look like:
  >transcript:ENSGXNT00000003375.p1 GENE.x~~y  ORF type:complete (+),psauron_score=0.964 len:104 ...
The protein name is the input FASTA record name with a trailing `.p<N>` ORF
index appended, so stripping `.p<N>` recovers the exact FASTA seqname. This
works for both renamed input FASTAs (`tid|type|species`) and bare-`tid`
projected FASTAs — no delimiter-specific parsing.

Outputs:
  --kept-fasta   FASTA of records NOT called coding (the retained non-coding set)
  --coding-ids   newline-separated list of coding seqnames (dropped)
  --mqc          one-row MultiQC custom-content TSV (n_in/n_coding/n_kept/pct)
"""
import argparse
import re
import sys
from typing import List, Set

PEP_HEADER_RE = re.compile(r"^>(\S+)")
ORF_TYPE_RE = re.compile(r"ORF\s+type[:=]\s*(\S+)", re.IGNORECASE)
PSAURON_RE = re.compile(r"psauron_score=([0-9.eE+-]+)")
# Strip a trailing ORF index ".p1", ".p2", ... to recover the FASTA seqname.
ORF_SUFFIX_RE = re.compile(r"\.p\d+$")


def parse_fasta_names(path: str) -> List[str]:
    """Return the ordered list of record names (token after '>')."""
    names: List[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
    return names


def iter_fasta_records(path: str):
    """Yield (name, header_line, [sequence_lines]) for each FASTA record."""
    name = None
    header = None
    seq: List[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, header, seq
                header = line
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, header, seq


def parse_coding(pep_path: str, psauron_min: float) -> Set[str]:
    """Return the set of FASTA seqnames with a complete, PSAURON-confident ORF."""
    coding: Set[str] = set()
    with open(pep_path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = PEP_HEADER_RE.match(line)
            if not m:
                continue
            seqname = ORF_SUFFIX_RE.sub("", m.group(1))

            orf_m = ORF_TYPE_RE.search(line)
            orf_type = orf_m.group(1).lower() if orf_m else ""
            is_complete = "complete" in orf_type

            ps_m = PSAURON_RE.search(line)
            try:
                psauron = float(ps_m.group(1)) if ps_m else 0.0
            except ValueError:
                psauron = 0.0

            if is_complete and psauron >= psauron_min:
                coding.add(seqname)
    return coding


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in-fasta", required=True, help="input lncRNA FASTA (TD2 was run on this)")
    p.add_argument("--pep", required=True, help="TD2.Predict .pep file")
    p.add_argument("--psauron-min", type=float, default=0.5,
                   help="minimum psauron_score for a complete ORF to count as coding (default 0.5)")
    p.add_argument("--kept-fasta", required=True, help="output FASTA of retained (non-coding) records")
    p.add_argument("--coding-ids", required=True, help="output list of dropped (coding) seqnames")
    p.add_argument("--mqc", required=True, help="output MultiQC custom-content TSV (one row)")
    p.add_argument("--sample", required=True, help="sample name for the MultiQC row")
    p.add_argument("--stage", required=True, choices=["input", "projected"],
                   help="pipeline stage (for the MultiQC table)")
    p.add_argument("--feature", default="lncRNA", help="feature label for the MultiQC row")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    names = parse_fasta_names(args.in_fasta)
    n_in = len(names)
    coding = parse_coding(args.pep, args.psauron_min)
    # Only count coding hits that actually correspond to an input record.
    name_set = set(names)
    coding = {c for c in coding if c in name_set}
    n_coding = len(coding)
    n_kept = n_in - n_coding
    pct = (100.0 * n_coding / n_in) if n_in else 0.0

    # Write kept (non-coding) FASTA, preserving input order and full headers.
    # FASTA filtering keys on the FULL record name (exact match).
    with open(args.kept_fasta, "w", encoding="utf-8") as out:
        for name, header, seq in iter_fasta_records(args.in_fasta):
            if name in coding:
                continue
            if header is not None:
                out.write(header)
            out.writelines(seq)

    # GFF subsetting downstream keys on the BARE transcript id, so emit the
    # token before the first '|' (renamed input headers are tid|type|species;
    # projected headers are already bare tid → split is a no-op there).
    bare_ids = sorted({c.split("|", 1)[0] for c in coding})
    with open(args.coding_ids, "w", encoding="utf-8") as out:
        for c in bare_ids:
            out.write(c + "\n")

    # MultiQC custom-content single-row table. Routing is by the *_td2_coding_mqc.tsv
    # sp pattern (defined in assets/multiqc/sections.yml); the file deliberately
    # carries NO `# id:` header — embedding one would double-render the section.
    with open(args.mqc, "w", encoding="utf-8") as out:
        out.write("Sample\tstage\tfeature\tn_in\tn_coding\tn_kept\tpct_coding\n")
        out.write(f"{args.sample}\t{args.stage}\t{args.feature}\t{n_in}\t{n_coding}\t{n_kept}\t{pct:.2f}\n")

    print(f"td2_coding_filter: {args.sample} ({args.stage}) "
          f"n_in={n_in} coding={n_coding} kept={n_kept} ({pct:.1f}% coding)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
