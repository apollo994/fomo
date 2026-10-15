#!/usr/bin/env python3
"""
Filter coding-potential transcripts out of a lncRNA FASTA using TD2 output.

A transcript is coding (and dropped) when it has at least one ORF in the TD2.Predict `.pep` file
that is both complete (ORF type contains "complete") and PSAURON-confident
(psauron_score >= --psauron-min).

TD2.Predict `.pep` headers look like:
  >transcript:ENSGXNT00000003375.p1 GENE.x~~y  ORF type:complete (+),psauron_score=0.964 len:104 ...
The protein name is the FASTA record name plus a `.p<N>` ORF index, so stripping `.p<N>` recovers
the record name exactly, whatever its delimiters (`tid|gene_class|species` for the input stage,
`tid|source` for the projected stage).

Outputs (--prefix P):
  P.noncoding.fasta      the records not called coding, in input order
  P.coding_ids.txt       the coding record names: the bare transcript id for --stage input
                         (informational), the name verbatim for --stage projected (FILTER_ALLMODELS
                         drops by the source-qualified id; a bare id could name another source's
                         transcript)
  P.td2_coding.data.json n_in / n_coding / n_kept as a `td2_coding` data record (with --identity;
                         the record's identity gets stage = --stage)
"""

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fomo_io  # noqa: E402
import fomo_stats  # noqa: E402

PEP_HEADER_RE = re.compile(r'^>(\S+)')
ORF_TYPE_RE = re.compile(r'ORF\s+type[:=]\s*(\S+)', re.IGNORECASE)
PSAURON_RE = re.compile(r'psauron_score=([^\s,;]*)')
ORF_SUFFIX_RE = re.compile(r'\.p\d+$')
log = fomo_io.get_logger('td2_coding_filter')


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--fasta', required=True, help='lncRNA FASTA that TD2 ran on')
    p.add_argument('--pep', required=True, help='TD2.Predict .pep file')
    p.add_argument('--psauron-min', type=float, default=0.5, help='minimum psauron_score of a coding ORF (0.5)')
    p.add_argument('--prefix', required=True, help='output path prefix')
    p.add_argument('--stage', required=True, choices=['input', 'projected'], help='pipeline stage')
    p.add_argument('--identity', help='identity of the data record (JSON object); without it no record')
    return p.parse_args()


def psauron_score(header: str) -> float | None:
    """The header's psauron_score, or None if it has no parseable one."""
    match = PSAURON_RE.search(header)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def coding_names(pep: str, psauron_min: float) -> set[str]:
    """Record names with a complete ORF whose psauron_score is >= psauron_min. A complete ORF
    without a parseable score is not coding, with a warning."""
    coding: set[str] = set()
    with fomo_io.open_text(pep) as fh:
        for line in fh:
            header = PEP_HEADER_RE.match(line)
            orf = ORF_TYPE_RE.search(line)
            if not header or not orf or 'complete' not in orf.group(1).lower():
                continue
            name = ORF_SUFFIX_RE.sub('', header.group(1))
            score = psauron_score(line)
            if score is None:
                log.warning(f'{pep}: no parseable psauron_score for {name}; not counted as coding')
            elif score >= psauron_min:
                coding.add(name)
    return coding


def write_noncoding(fasta: str, coding: set[str], output: str) -> tuple[int, set[str]]:
    """Copy the non-coding records; -> (records read, coding names found in the FASTA)."""
    n_in, found = 0, set()
    with fomo_io.open_text(output, 'wt') as out:
        for name, header, seq in fomo_io.iter_fasta(fasta):
            n_in += 1
            if name in coding:
                found.add(name)
                continue
            out.write(header)
            out.writelines(seq)
    return n_in, found


def main() -> int:
    args = parse_args()
    coding = coding_names(args.pep, args.psauron_min)
    n_in, found = write_noncoding(args.fasta, coding, f'{args.prefix}.noncoding.fasta')
    n_coding = len(found)
    n_kept = n_in - n_coding
    pct = 100.0 * n_coding / n_in if n_in else 0.0

    drop_ids = sorted({c.split('|', 1)[0] for c in found}) if args.stage == 'input' else sorted(found)
    with fomo_io.open_text(f'{args.prefix}.coding_ids.txt', 'wt') as out:
        out.writelines(c + '\n' for c in drop_ids)

    if args.identity:
        identity = fomo_stats.identity_arg(args.identity) | {'stage': args.stage}
        data = {'n_in': n_in, 'n_coding': n_coding, 'n_kept': n_kept, 'psauron_min': args.psauron_min}
        fomo_stats.write_records(
            f'{args.prefix}.td2_coding.data.json', fomo_stats.make_record('td2_coding', identity, data)
        )

    log.info(f'{args.prefix} ({args.stage}) n_in={n_in} coding={n_coding} kept={n_kept} ({pct:.1f}% coding)')
    return 0


if __name__ == '__main__':
    fomo_io.run(main)
