#!/usr/bin/env python3
"""
Convert AGAT `agat_sp_statistics.pl --yaml` output into a single MultiQC
custom-content TSV row.

Curated columns:
  role, n_genes, n_transcripts, n_exons, n_introns,
  mean_transcript_length, min_transcript_length, max_transcript_length,
  mean_exon_length, max_exon_length,
  mean_intron_length, min_intron_length, max_intron_length,
  genome_coverage

The section ID is set via --section-id (e.g. `agat_input` or
`agat_projection`) so the same script feeds multiple MultiQC sections.

AGAT YAML keys are case-insensitive and may appear with/without an `s`
suffix on the section name (`with_isoforms` / `without_isoform` /
`without_isoforms`). This script:

  1. Loads the YAML
  2. Flattens it into `section__metric` string keys (lowercased,
     punctuation stripped) so regex matching is robust
  3. Aggregates across feature types (sums for counts, weighted means
     for means, min/max across features for extrema)
  4. Resolves `genome_coverage` from AGAT's `% of the genome covered by
     <tag>` field (only present when --gs was passed to AGAT)
"""
import argparse
import re
import sys
from typing import Any, Dict, List, Optional

import yaml


def normalise_key(s: str) -> str:
    """Lowercase + whitespace→underscore + strip non-[a-z0-9_]."""
    s = str(s).strip().lower()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^a-z0-9_]", "", s)
    return s


def flatten(node: Any, prefix: str = "") -> Dict[str, Any]:
    """Recursively flatten a nested dict into 'a__b__c' keys."""
    out: Dict[str, Any] = {}
    if isinstance(node, dict):
        for k, v in node.items():
            key = normalise_key(k)
            new_prefix = f"{prefix}__{key}" if prefix else key
            out.update(flatten(v, new_prefix))
    elif isinstance(node, list):
        for i, v in enumerate(node):
            out.update(flatten(v, f"{prefix}__{i}"))
    else:
        out[prefix] = node
    return out


def find_matches(flat: Dict[str, Any], patterns: List[str]) -> List[Any]:
    """Return values whose keys match the first regex that finds any."""
    for pat in patterns:
        regex = re.compile(pat, re.IGNORECASE)
        values = [v for k, v in flat.items() if regex.search(k)]
        if values:
            return values
    return []


def agg_sum(values: List[Any]) -> Optional[float]:
    nums = [v for v in values if isinstance(v, (int, float)) and not isinstance(v, bool)]
    return sum(nums) if nums else None


def agg_min(values: List[Any]) -> Optional[float]:
    nums = [v for v in values if isinstance(v, (int, float)) and not isinstance(v, bool)]
    return min(nums) if nums else None


def agg_max(values: List[Any]) -> Optional[float]:
    nums = [v for v in values if isinstance(v, (int, float)) and not isinstance(v, bool)]
    return max(nums) if nums else None


def fmt(v: Any) -> str:
    if v is None:
        return "NA"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        if v.is_integer():
            return str(int(v))
        return f"{v:.2f}"
    return str(v)


# Patterns probed in order; first regex with any match wins.
PATTERNS: Dict[str, List[str]] = {
    "n_genes":               [r"__value__number_of_(ncrna_)?gene\b"],
    "n_transcripts":         [r"__value__number_of_(mrna|lnc_rna|transcript)\b"],
    "n_exons":               [r"__value__number_of_exon\b"],
    "n_introns":             [r"__value__number_of_intron_in_exon\b"],
    "max_transcript_length": [r"__value__longest_(mrna|lnc_rna|transcript)_bp"],
    "min_transcript_length": [r"__value__shortest_(mrna|lnc_rna|transcript)_bp"],
    "max_exon_length":       [r"__value__longest_exon_bp"],
    "max_intron_length":     [r"__value__longest_intron_into_exon_part_bp"],
    "min_intron_length":     [r"__value__shortest_intron_into_exon_part_bp"],
}


def compute_metric(flat: Dict[str, Any], column: str) -> Optional[float]:
    """Aggregate the raw matches for `column` into a final scalar."""
    if column.startswith("n_"):
        return agg_sum(find_matches(flat, PATTERNS[column]))
    if column.startswith("max_"):
        return agg_max(find_matches(flat, PATTERNS[column]))
    if column.startswith("min_"):
        return agg_min(find_matches(flat, PATTERNS[column]))

    if column == "mean_transcript_length":
        total = agg_sum(find_matches(flat, [r"__value__total_(mrna|lnc_rna|transcript)_length_bp"]))
        n     = agg_sum(find_matches(flat, [r"__value__number_of_(mrna|lnc_rna|transcript)\b"]))
        return (total / n) if (total is not None and n) else None

    if column == "mean_exon_length":
        total = agg_sum(find_matches(flat, [r"__value__total_exon_length_bp"]))
        n     = agg_sum(find_matches(flat, [r"__value__number_of_exon\b"]))
        return (total / n) if (total is not None and n) else None

    if column == "mean_intron_length":
        total = agg_sum(find_matches(flat, [r"__value__total_intron_length_per_exon_bp"]))
        n     = agg_sum(find_matches(flat, [r"__value__number_of_intron_in_exon\b"]))
        return (total / n) if (total is not None and n) else None

    if column == "genome_coverage":
        # AGAT (--gs) emits `% of genome covered by <tag>` keys. After
        # normalisation `%` is stripped, leaving `_of_genome_covered_by_...`
        # glued onto `value` by `__` → 3+ underscores between them. Prefer
        # the broadest feature type when several are present.
        for tag in ("gene", "ncrna_gene", "mrna", "lnc_rna", "transcript"):
            matches = find_matches(flat, [rf"__value_+of_genome_covered_by_{tag}\b"])
            if matches:
                return agg_max(matches)
        return None

    return None


COLUMNS: List[str] = [
    "role",
    "n_genes",
    "n_transcripts",
    "n_exons",
    "n_introns",
    "mean_transcript_length",
    "min_transcript_length",
    "max_transcript_length",
    "mean_exon_length",
    "max_exon_length",
    "mean_intron_length",
    "min_intron_length",
    "max_intron_length",
    "genome_coverage",
]


SECTION_NAMES: Dict[str, str] = {
    "agat_input":      "AGAT statistics — input annotations",
    "agat_projection": "AGAT statistics — projected models",
}


def write_tsv(out_path: str, sample_name: str, role: str,
              section_id: str, flat: Dict[str, Any]) -> int:
    """Write the one-row TSV. Returns count of resolved (non-NA) numeric columns.

    No `# id:` header block is written — section routing is handled by the
    `sp.<section_id>.fn` pattern in MultiQC config. Embedding a `# id:`
    header here would cause MultiQC to render the section twice (once via
    the `sp` match, once via the embedded id)."""
    _ = section_id  # accepted for symmetry with module signature; not used
    resolved: Dict[str, str] = {"role": role}
    for col in COLUMNS:
        if col == "role":
            continue
        resolved[col] = fmt(compute_metric(flat, col))

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("Sample\t" + "\t".join(COLUMNS) + "\n")
        fh.write(sample_name + "\t" + "\t".join(resolved[c] for c in COLUMNS) + "\n")

    return sum(1 for c, v in resolved.items() if c != "role" and v != "NA")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",      required=True, metavar="YAML",
                   help="AGAT YAML stats output")
    p.add_argument("--name",       required=True, metavar="SAMPLE",
                   help="Row sample identifier")
    p.add_argument("--role",       required=True, choices=["source", "target", "projection"],
                   help="Role label for the row")
    p.add_argument("--section-id", required=True, choices=list(SECTION_NAMES),
                   help="MultiQC section id (agat_input or agat_projection)")
    p.add_argument("--output",     required=True, metavar="TSV",
                   help="Output one-row TSV")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    with open(args.input, "r", encoding="utf-8") as fh:
        try:
            data = yaml.safe_load(fh)
        except yaml.YAMLError as e:
            print(f"ERROR: failed to parse YAML {args.input}: {e}", file=sys.stderr)
            return 1
    if data is None:
        data = {}

    flat = flatten(data)
    n_resolved = write_tsv(args.output, args.name, args.role, args.section_id, flat)
    print(
        f"AGAT_TO_MQC: {args.name} role={args.role} section={args.section_id} → "
        f"{n_resolved}/{len(COLUMNS) - 1} numeric metrics resolved",
        file=sys.stderr,
    )
    if n_resolved == 0:
        sample_keys = sorted(flat.keys())[:10]
        print(
            f"WARNING: no metrics extracted. First flat keys: {sample_keys}",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
