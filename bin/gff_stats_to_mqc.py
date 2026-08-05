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
from typing import Any, Dict, List, Optional

NA = "NA"

# Column order for the per-transcript-type table.
TRANSCRIPT_COLUMNS: List[str] = [
    "role",
    "kind",
    "transcript_type",
    "n_genes",
    "n_transcripts",
    "n_exons",
    "mean_exon_length",
    "max_exon_length",
    "n_introns",
    "mean_intron_length",
    "min_intron_length",
    "max_intron_length",
    "mean_transcript_length",
    "min_transcript_length",
    "max_transcript_length",
    "mean_spliced_length",
    "min_spliced_length",
    "max_spliced_length",
    "n_cds",
    "mean_cds_length",
    "top_biotype",
]

# Column order for the per-gene-category table.
GENE_COLUMNS: List[str] = [
    "role",
    "kind",
    "gene_category",
    "n_genes",
    "mean_gene_length",
    "min_gene_length",
    "max_gene_length",
    "n_transcript_types",
    "top_biotype",
]


def fmt(v: Any) -> str:
    """Render a metric for the TSV. None → NA; integral floats lose the .0."""
    if v is None:
        return NA
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return str(int(v)) if v.is_integer() else f"{v:.2f}"
    return str(v)


def length_field(block: Optional[Dict[str, Any]], key: str) -> Optional[Any]:
    """Pull min/max/mean out of a `{min, max, mean}` sub-object, or None."""
    if not isinstance(block, dict):
        return None
    value = block.get(key)
    return value if isinstance(value, (int, float)) else None


def top_key(counts: Any) -> Optional[str]:
    """Highest-count key of a `{name: count}` map, ties broken by name."""
    if not isinstance(counts, dict) or not counts:
        return None
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def transcript_rows(stats: Dict[str, Any], role: str, kind: str) -> List[Dict[str, str]]:
    """Build one row per transcript type, in the JSON's own (count-desc) order."""
    rows: List[Dict[str, str]] = []
    for ttype, ts in (stats.get("transcript_type_stats") or {}).items():
        length_stats = ts.get("length_stats")
        exon = ts.get("exon_stats") or {}
        exon_len = exon.get("length")
        spliced = exon.get("concatenated_length")
        intron = ts.get("intron_stats") or {}
        intron_len = intron.get("length")
        cds = ts.get("cds_stats") or {}
        cds_len = cds.get("length")
        assoc = ts.get("associated_genes") or {}

        rows.append(
            {
                "role": role,
                "kind": kind,
                "transcript_type": ttype,
                "n_genes": fmt(assoc.get("total_count")),
                "n_transcripts": fmt(ts.get("total_count")),
                "n_exons": fmt(exon.get("total_count")),
                "mean_exon_length": fmt(length_field(exon_len, "mean")),
                "max_exon_length": fmt(length_field(exon_len, "max")),
                "n_introns": fmt(intron.get("total_count")),
                "mean_intron_length": fmt(length_field(intron_len, "mean")),
                "min_intron_length": fmt(length_field(intron_len, "min")),
                "max_intron_length": fmt(length_field(intron_len, "max")),
                "mean_transcript_length": fmt(length_field(length_stats, "mean")),
                "min_transcript_length": fmt(length_field(length_stats, "min")),
                "max_transcript_length": fmt(length_field(length_stats, "max")),
                "mean_spliced_length": fmt(length_field(spliced, "mean")),
                "min_spliced_length": fmt(length_field(spliced, "min")),
                "max_spliced_length": fmt(length_field(spliced, "max")),
                "n_cds": fmt(cds.get("total_count")),
                "mean_cds_length": fmt(length_field(cds_len, "mean")),
                "top_biotype": fmt(top_key(ts.get("biotype_counts"))),
            }
        )

    if not rows:
        placeholder = {c: NA for c in TRANSCRIPT_COLUMNS}
        placeholder.update({"role": role, "kind": kind, "transcript_type": "none"})
        rows.append(placeholder)
    return rows


def gene_rows(stats: Dict[str, Any], role: str, kind: str) -> List[Dict[str, str]]:
    """Build one row per gene category (coding / pseudogene / non_coding)."""
    rows: List[Dict[str, str]] = []
    for category, gc in (stats.get("gene_category_stats") or {}).items():
        length_stats = gc.get("length_stats")
        ttype_counts = gc.get("transcript_type_counts")
        rows.append(
            {
                "role": role,
                "kind": kind,
                "gene_category": category,
                "n_genes": fmt(gc.get("total_count")),
                "mean_gene_length": fmt(length_field(length_stats, "mean")),
                "min_gene_length": fmt(length_field(length_stats, "min")),
                "max_gene_length": fmt(length_field(length_stats, "max")),
                "n_transcript_types": fmt(
                    len(ttype_counts) if isinstance(ttype_counts, dict) else None
                ),
                "top_biotype": fmt(top_key(gc.get("biotype_counts"))),
            }
        )

    if not rows:
        placeholder = {c: NA for c in GENE_COLUMNS}
        placeholder.update({"role": role, "kind": kind, "gene_category": "none"})
        rows.append(placeholder)
    return rows


def write_tsv(path: str, base_name: str, name_suffix_col: str,
              columns: List[str], rows: List[Dict[str, str]]) -> None:
    """Write the table. Sample name = `<base>.<value of name_suffix_col>`."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("Sample\t" + "\t".join(columns) + "\n")
        for row in rows:
            sample = f"{base_name}.{row[name_suffix_col]}"
            fh.write(sample + "\t" + "\t".join(row[c] for c in columns) + "\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, metavar="JSON",
                   help="gff-feature-stats JSON output")
    p.add_argument("--name", required=True, metavar="SAMPLE",
                   help="Row sample identifier; the transcript type / gene "
                        "category is appended to it")
    p.add_argument("--role", required=True,
                   choices=["source", "target", "projection"],
                   help="Role label for the rows")
    p.add_argument("--kind", default=NA,
                   help="Annotation stage (raw/filtered/decoy/projected)")
    p.add_argument("--section-id", required=True,
                   choices=["gffstats_input", "gffstats_projection"],
                   help="MultiQC section id of the transcript table")
    p.add_argument("--output", required=True, metavar="TSV",
                   help="Output transcript-type TSV")
    p.add_argument("--genes-output", metavar="TSV",
                   help="Optional output gene-category TSV")
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
        f"section={args.section_id} → {len(t_rows)} transcript-type row(s), "
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
