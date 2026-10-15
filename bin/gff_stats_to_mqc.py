#!/usr/bin/env python3
"""
Convert a `gff-feature-stats` JSON document into MultiQC custom-content TSVs.

Two tables are produced from one JSON:

  * transcript table (always) — one row per transcript type present in
    `transcript_type_stats`, sample name `<--name>.<transcript_type>`.
  * gene table (only with --genes-output) — one row per category in
    `gene_category_stats`, sample name `<--name>.<category>`.

Rows are emitted even when the JSON is empty (an all-NA placeholder row), so
MultiQC never has to parse a header-only file.

No `# id:` header block is written — section routing is handled entirely by the
`sp.<section>.fn` patterns in assets/multiqc/sections.yml. Embedding a `# id:`
here would make MultiQC discover the section twice and render it twice.

Blocks that gff-feature-stats omits are reported as NA rather than 0, because
absence is meaningful: `intron_stats` is missing for a type whose transcripts
are all single-exon, `exon_stats.concatenated_length` for a type with no
multi-exon transcript, and `cds_stats` for a non-coding type.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fomo_stats  # noqa: E402

NA = fomo_stats.NA
TRANSCRIPT_COLUMNS = fomo_stats.GFFSTATS_TRANSCRIPT_COLUMNS
GENE_COLUMNS = fomo_stats.GFFSTATS_GENE_COLUMNS
transcript_rows = fomo_stats.gffstats_transcript_rows
gene_rows = fomo_stats.gffstats_gene_rows


def write_tsv(path: str, base_name: str, name_suffix_col: str, columns: List[str], rows: List[Dict[str, str]]) -> None:
    """Write the table. Sample name = `<base>.<value of name_suffix_col>`."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("Sample\t" + "\t".join(columns) + "\n")
        for row in rows:
            sample = f"{base_name}.{row[name_suffix_col]}"
            fh.write(sample + "\t" + "\t".join(row[c] for c in columns) + "\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, metavar="JSON", help="gff-feature-stats JSON output")
    p.add_argument(
        "--name",
        required=True,
        metavar="SAMPLE",
        help="Row sample identifier; the transcript type / gene category is appended to it",
    )
    p.add_argument("--role", required=True, choices=["source", "target", "projection"], help="Role label for the rows")
    p.add_argument("--kind", default=NA, help="Annotation stage (raw/filtered/decoy/projected)")
    p.add_argument("--output", required=True, metavar="TSV", help="Output transcript-type TSV")
    p.add_argument("--genes-output", metavar="TSV", help="Optional output gene-category TSV")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    with open(args.input, "r", encoding="utf-8") as fh:
        try:
            stats = json.load(fh)
        except json.JSONDecodeError as e:
            print(f"ERROR: failed to parse JSON {args.input}: {e}", file=sys.stderr)
            return 1
    if not isinstance(stats, dict):
        stats = {}

    t_rows = transcript_rows(stats, args.role, args.kind)
    write_tsv(args.output, args.name, "transcript_type", TRANSCRIPT_COLUMNS, t_rows)

    n_na = sum(1 for r in t_rows for c, v in r.items() if c in TRANSCRIPT_COLUMNS and v == NA)
    n_cells = len(t_rows) * (len(TRANSCRIPT_COLUMNS) - 3)  # minus the 3 label columns
    print(
        f"GFF_STATS_TO_MQC: {args.name} role={args.role} kind={args.kind} "
        f"→ {len(t_rows)} transcript-type row(s), "
        f"{n_cells - n_na}/{n_cells} metrics resolved",
        file=sys.stderr,
    )

    if args.genes_output:
        g_rows = gene_rows(stats, args.role, args.kind)
        write_tsv(args.genes_output, args.name, "gene_category", GENE_COLUMNS, g_rows)
        print(
            f"GFF_STATS_TO_MQC: {args.name} → {len(g_rows)} gene-category row(s)",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
