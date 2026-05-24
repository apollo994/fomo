#!/usr/bin/env python3
"""
Convert AGAT `agat_sp_statistics.pl --yaml` output into MultiQC
custom-content TSV files.

Produces TWO TSVs per input YAML:

  1. <prefix>_summary_mqc.tsv  - curated 8-metric headline table; one row
                                 per AGAT sample, MultiQC id `agat_summary`.
  2. <prefix>_full_mqc.tsv     - exhaustive per-feature-type table; one row
                                 per top-level feature found in the YAML
                                 (e.g. `lnc_rna`, `mrna`), with every
                                 numeric metric AGAT emits as a column.
                                 MultiQC id `agat_full`.

The AGAT YAML schema varies — keys may appear in any case, with or without
trailing 's', and the section name alternates between `with_isoforms` /
`without_isoform` / `without_isoforms`. This parser:

  1. Loads the YAML
  2. Flattens it into `section__metric` string keys (lowercased,
     punctuation stripped) so regex matching is robust
  3. For the summary: matches a curated list of metrics by regex,
     aggregating sensibly across feature types
  4. For the full table: emits one row per top-level feature, with every
     numeric metric (path-stripped) as a column

Missing curated metrics are written as "NA". A summary line is printed
to stderr.
"""
import argparse
import re
import sys
from typing import Any, Dict, List, Optional, Tuple

import yaml


# ---------------------------------------------------------------------------
# Curated summary metrics: (column header, list of regex patterns).
# Patterns try in order; first pattern that matches at least one flat key
# wins (so `with_isoforms` is preferred over `without_isoforms?`).
# ---------------------------------------------------------------------------
METRIC_SPECS: List[Tuple[str, List[str]]] = [
    ("n_genes",                 [r"with_isoforms__value__number_of_(ncrna_)?gene\b",
                                 r"without_isoforms?__value__number_of_(ncrna_)?gene\b",
                                 r"__value__number_of_(ncrna_)?gene\b"]),
    ("n_transcripts",           [r"with_isoforms__value__number_of_(mrna|lnc_rna|transcript)\b",
                                 r"without_isoforms?__value__number_of_(mrna|lnc_rna|transcript)\b",
                                 r"__value__number_of_(mrna|lnc_rna|transcript)\b"]),
    ("n_exons",                 [r"with_isoforms__value__number_of_exon\b",
                                 r"without_isoforms?__value__number_of_exon\b",
                                 r"__value__number_of_exon\b"]),
    ("n_introns",               [r"with_isoforms__value__number_of_intron_in_exon",
                                 r"without_isoforms?__value__number_of_intron_in_exon"]),
    ("n_cds",                   [r"with_isoforms__value__number_of_cds\b",
                                 r"without_isoforms?__value__number_of_cds\b"]),
    ("mean_transcript_length",  [r"with_isoforms__value__mean_(mrna|lnc_rna)_length",
                                 r"without_isoforms?__value__mean_(mrna|lnc_rna)_length"]),
    ("total_transcript_length", [r"with_isoforms__value__total_(mrna|lnc_rna)_length",
                                 r"without_isoforms?__value__total_(mrna|lnc_rna)_length"]),
    ("mean_exon_length",        [r"with_isoforms__value__mean_exon_length",
                                 r"without_isoforms?__value__mean_exon_length"]),
]

# Inside-feature path that hosts the actual metrics. Tried in order; first
# match wins. Covers all three AGAT layout variants.
FEATURE_VALUE_PATHS = ["with_isoforms__value", "without_isoforms__value", "without_isoform__value", "value"]


def normalise_key(s: str) -> str:
    """Lowercase + whitespace→underscore + strip non-[a-z0-9_] characters."""
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


def find_all_matches(flat: Dict[str, Any], patterns: List[str]) -> List[Any]:
    """Return values whose keys match the first pattern that finds any."""
    for pat in patterns:
        regex = re.compile(pat, re.IGNORECASE)
        values = [v for k, v in flat.items() if regex.search(k)]
        if values:
            return values
    return []


def aggregate_sum(values: List[Any]) -> Optional[float]:
    if not values:
        return None
    total = 0.0
    for v in values:
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            total += v
        else:
            return None
    return total


def format_value(v: Any) -> str:
    if v is None:
        return "NA"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        if v.is_integer():
            return str(int(v))
        return f"{v:.2f}"
    return str(v)


# ---------------------------------------------------------------------------
# Summary table writer (curated metrics, single row per sample)
# ---------------------------------------------------------------------------
def write_summary_tsv(out_path: str, sample_name: str, flat: Dict[str, Any]) -> int:
    """Write the curated 8-metric summary TSV. Returns count of resolved
    (non-NA) metrics."""
    TOTAL_TRANSCRIPT_LEN_PATTERNS = [
        r"with_isoforms__value__total_(mrna|lnc_rna)_length",
        r"without_isoforms?__value__total_(mrna|lnc_rna)_length",
    ]
    TOTAL_EXON_LEN_PATTERNS = [
        r"with_isoforms__value__total_exon_length",
        r"without_isoforms?__value__total_exon_length",
    ]

    resolved: Dict[str, str] = {}
    for column, patterns in METRIC_SPECS:
        matches = find_all_matches(flat, patterns)
        if column == "mean_transcript_length":
            n = aggregate_sum(find_all_matches(flat, [
                r"with_isoforms__value__number_of_(mrna|lnc_rna|transcript)\b",
                r"without_isoforms?__value__number_of_(mrna|lnc_rna|transcript)\b",
            ]))
            total = aggregate_sum(find_all_matches(flat, TOTAL_TRANSCRIPT_LEN_PATTERNS))
            value = (total / n) if (n and total is not None) else None
        elif column == "mean_exon_length":
            n = aggregate_sum(find_all_matches(flat, [
                r"with_isoforms__value__number_of_exon\b",
                r"without_isoforms?__value__number_of_exon\b",
            ]))
            total = aggregate_sum(find_all_matches(flat, TOTAL_EXON_LEN_PATTERNS))
            value = (total / n) if (n and total is not None) else None
        else:
            value = aggregate_sum(matches)
        resolved[column] = format_value(value)

    columns = [c for c, _ in METRIC_SPECS]
    header_block = [
        "# id: agat_summary",
        "# section_name: 'AGAT annotation statistics (summary)'",
        "# description: 'Headline counts and aggregate lengths from agat_sp_statistics --yaml. See AGAT full table for per-feature-type detail.'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'agat_summary_table'",
        "#     namespace: 'AGAT'",
    ]
    with open(out_path, "w", encoding="utf-8") as fh:
        for line in header_block:
            fh.write(line + "\n")
        fh.write("Sample\t" + "\t".join(columns) + "\n")
        fh.write(sample_name + "\t" + "\t".join(resolved[c] for c in columns) + "\n")

    return sum(1 for v in resolved.values() if v != "NA")


# ---------------------------------------------------------------------------
# Full table writer (per-feature-type rows, every numeric metric)
# ---------------------------------------------------------------------------
def extract_feature_metrics(data: Any) -> Dict[str, Dict[str, Any]]:
    """For each top-level feature (e.g. `mrna`, `lnc_rna`), return a dict of
    {metric_name: value} extracted from the deepest available value path
    (`with_isoforms.value` preferred, then `without_isoforms.value` /
    `without_isoform.value` / `value`)."""
    if not isinstance(data, dict):
        return {}
    out: Dict[str, Dict[str, Any]] = {}
    for feat_name, feat_data in data.items():
        if not isinstance(feat_data, dict):
            continue
        # Find the value subtree
        value_dict: Optional[Dict[str, Any]] = None
        for candidate in ["with_isoforms", "without_isoforms", "without_isoform"]:
            sub = feat_data.get(candidate)
            if isinstance(sub, dict) and isinstance(sub.get("value"), dict):
                value_dict = sub["value"]
                break
        if value_dict is None:
            # Fallback: direct `value` at feature level
            if isinstance(feat_data.get("value"), dict):
                value_dict = feat_data["value"]
        if value_dict is None:
            continue

        metrics: Dict[str, Any] = {}
        for k, v in value_dict.items():
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                metrics[normalise_key(k)] = v
        if metrics:
            out[normalise_key(feat_name)] = metrics
    return out


def write_full_tsv(out_path: str, sample_name: str, data: Any) -> int:
    """Write the per-feature-type full metrics TSV. Returns row count."""
    feature_metrics = extract_feature_metrics(data)
    if not feature_metrics:
        # Empty: still write a stub so MultiQC has something
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write("# id: agat_full\n")
            fh.write("# section_name: 'AGAT annotation statistics (full)'\n")
            fh.write("# plot_type: 'table'\n")
            fh.write("Sample\tnote\n")
            fh.write(f"{sample_name}\tno_metrics_found\n")
        return 0

    # Union of all metric keys across feature types — gives the column set
    all_metrics: List[str] = sorted({m for feat in feature_metrics.values() for m in feat.keys()})

    header_block = [
        "# id: agat_full",
        "# section_name: 'AGAT annotation statistics (full)'",
        "# description: 'Exhaustive per-feature-type metrics from agat_sp_statistics --yaml. One row per top-level feature type (mrna, lnc_rna, gene, ...). Columns flatten AGAT metric names. Click table column headers to sort.'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'agat_full_table'",
        "#     namespace: 'AGAT'",
        "#     only_defined_headers: false",
    ]

    with open(out_path, "w", encoding="utf-8") as fh:
        for line in header_block:
            fh.write(line + "\n")
        fh.write("Sample\tfeature\t" + "\t".join(all_metrics) + "\n")
        for feat, metrics in sorted(feature_metrics.items()):
            row_id = f"{sample_name}__{feat}"
            cells = [format_value(metrics.get(m)) for m in all_metrics]
            fh.write(row_id + "\t" + feat + "\t" + "\t".join(cells) + "\n")

    return len(feature_metrics)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",          required=True, metavar="YAML",
                   help="AGAT YAML stats output")
    p.add_argument("--name",           required=True, metavar="SAMPLE",
                   help="Sample name for the MultiQC row")
    p.add_argument("--summary-output", required=True, metavar="TSV",
                   help="Output curated summary TSV (one row, MultiQC id agat_summary)")
    p.add_argument("--full-output",    required=True, metavar="TSV",
                   help="Output full per-feature-type TSV (MultiQC id agat_full)")
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
        print(f"WARNING: {args.input} is empty or null", file=sys.stderr)
        data = {}

    flat = flatten(data)

    n_summary = write_summary_tsv(args.summary_output, args.name, flat)
    n_full    = write_full_tsv(args.full_output, args.name, data)

    print(
        f"AGAT_TO_MQC: {args.name} → summary={n_summary}/{len(METRIC_SPECS)} metrics, "
        f"full={n_full} feature-type rows",
        file=sys.stderr,
    )
    if n_summary == 0 and n_full == 0:
        print(
            f"WARNING: no metrics extracted. Top-level keys: {sorted(flat.keys())[:10]}",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
