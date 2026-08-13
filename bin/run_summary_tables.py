#!/usr/bin/env python3
"""
Digest a whole FOMO run into run-level MultiQC sections.

The per-target reports answer "how did this target do". This answers "what did this
run do": species per role, how much annotation went in and how much survived
filtering, how much each target gained, and which donors are worth using.

Everything is derived from artefacts the pipeline already produces — the
`*_mqc.tsv` adapter tables, the gffcompare `.stats` files and the top_sources CSVs
— so this script adds no new measurement, only cross-target arithmetic. Identity
comes from the `Sample` column (the ext.sample_name ladders in conf/modules.config
encode target/source/class) and from the stats filename grammar in fomo_stats.py.
Nothing here parses a GFF or a BAM.

Outputs, in report order (see assets/multiqc/run_summary_main.yml):
  run_overview_mqc.tsv            run overview: roles, pair counts, params
  run_species_mqc.tsv             per species: role + its own annotation input
  run_source_funnel_mqc.{tsv,json} per source: raw → spliced → non-coding funnel
  run_target_types_mqc.tsv        per (target × transcript type): before/added/after
  run_target_before_after_mqc.json  the same, as a stacked bargraph
  run_accuracy_heatmap_{f1,sn,pr}_mqc.json  source × target accuracy matrices
  run_aggregate_accuracy_mqc.{tsv,json}  aggregate models vs best single donor
  run_donor_ranking_mqc.tsv       per source: F1 across targets, top-N picks
  run_donor_picks_mqc.json        top-N picks per species, as a bargraph
  run_top_sources_mqc.tsv         the actual rank_1..rank_N per target
  run_summary.json                every number above, machine-readable
"""
import argparse
import csv
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence

# Nextflow only puts bin/ on PATH, not on Python's import path — but it stages the
# whole directory, so a __file__-relative insert finds the sibling module.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fomo_stats import AGGREGATE_IDS, f1, parse_identity, transcript_sn_pr  # noqa: E402

NA = "NA"

# (feature_type, decoy) → the `ext.sample_name` order token that projection and
# alignment rows carry. Defined in conf/modules.config, and deliberately sparse: a
# disabled track simply produces no row, so a missing token is normal, never an error.
# (The input-annotation side uses its own 1_raw / 2_lncRNA / 3_decoy_lncRNA /
# 4_mRNA / 5_decoy_mRNA ladder, read straight off the Sample column.)
PROJECTION_TOKEN = {
    ("lnc_RNA", False): "1_lncRNA",
    ("lnc_RNA", True):  "2_decoy_lncRNA",
    ("mRNA", False):    "3_mRNA",
    ("mRNA", True):     "4_decoy_mRNA",
}
# transcript_type → the real (non-decoy) token carrying its projected models.
TYPE_TO_TOKEN = {ft: PROJECTION_TOKEN[(ft, False)] for ft in ("lnc_RNA", "mRNA")}


# ── small helpers ───────────────────────────────────────────────────────────────
def num(value: Optional[str]) -> Optional[float]:
    """A metric cell as a float, or None for NA/missing/unparseable."""
    if value is None or value in (NA, "", "None"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def fmt(value) -> str:
    """Render a metric for a MultiQC TSV: None → NA, integral floats without .0."""
    if value is None:
        return NA
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        if value != value:                       # NaN
            return NA
        return str(int(value)) if value.is_integer() else f"{value:.2f}"
    return str(value)


def pct(part: Optional[float], whole: Optional[float]) -> Optional[float]:
    if part is None or whole in (None, 0):
        return None
    return 100.0 * part / whole


def add(*values: Optional[float]) -> Optional[float]:
    """Sum that propagates None — 'unknown + 5' is unknown, not 5."""
    if any(v is None for v in values):
        return None
    return sum(values)  # type: ignore[arg-type]


def split_sample(sample: str) -> Optional[tuple]:
    """`<prefix>.<order_token>.<type>` → (prefix, order_token, type).

    Species ids and transcript types contain no '.', so splitting from the right is
    exact. Returns None if the cell has too few parts to be one of ours.
    """
    parts = sample.split(".")
    if len(parts) < 3:
        return None
    return ".".join(parts[:-2]), parts[-2], parts[-1]


def read_tsvs(paths: Sequence[Path]) -> List[dict]:
    rows: List[dict] = []
    for path in paths:
        with path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                row["__file"] = path.name
                rows.append(row)
    return rows


def write_tsv(name: str, header: Sequence[str], rows: Sequence[Sequence]) -> None:
    with open(name, "w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(fmt(cell) for cell in row) + "\n")
    print(f"run_summary: wrote {name} ({len(rows)} rows)", file=sys.stderr)


def write_json(name: str, doc: dict) -> None:
    with open(name, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2)
    print(f"run_summary: wrote {name}", file=sys.stderr)


# ── input collection ────────────────────────────────────────────────────────────
class Inputs:
    """Everything the run produced, indexed by identity rather than by filename."""

    def __init__(self, mqc_dir: Path, roles_csv: Path, top_source_files: Sequence[Path]):
        files = sorted(p for p in mqc_dir.rglob("*") if p.is_file())
        self._assert_unique(files)

        def pick(suffix: str) -> List[Path]:
            return [p for p in files if p.name.endswith(suffix)]

        # role/kind/transcript-type tables of the INPUT annotations (source + target)
        self.input_rows = read_tsvs(pick("_gffstats_input_mqc.tsv"))
        self.gene_rows = read_tsvs(pick("_gffstats_genes_mqc.tsv"))
        self.projection_rows = read_tsvs(pick("_gffstats_projection_mqc.tsv"))
        self.td2_rows = read_tsvs(pick("_td2_coding_mqc.tsv"))
        self.align_rows = read_tsvs(pick("_samtools_align_mqc.tsv"))
        self.stats_files = pick(".gffcompare.stats")

        self.roles = self._read_roles(roles_csv)
        self.top_sources = self._read_top_sources(top_source_files)

        print(
            "run_summary: staged "
            f"{len(self.input_rows)} input, {len(self.gene_rows)} gene, "
            f"{len(self.projection_rows)} projection, {len(self.td2_rows)} TD2, "
            f"{len(self.align_rows)} alignment rows; "
            f"{len(self.stats_files)} gffcompare stats; "
            f"{len(self.top_sources)} top_sources",
            file=sys.stderr,
        )

    @staticmethod
    def _assert_unique(files: Sequence[Path]) -> None:
        """Basenames are unique across the whole pipeline by construction (every
        target-scoped ext.prefix starts with the target id, source-side ones are
        species-keyed). A duplicate means a naming contract broke upstream and two
        different measurements are about to be read as one — fail, loudly."""
        seen: Dict[str, Path] = {}
        for path in files:
            if path.name in seen:
                sys.exit(
                    f"ERROR: duplicate staged basename '{path.name}' "
                    f"({seen[path.name]} vs {path}). Two stat files claim the same "
                    "identity — check the ext.prefix closures in conf/modules.config."
                )
            seen[path.name] = path

    @staticmethod
    def _read_roles(path: Path) -> Dict[str, dict]:
        roles: Dict[str, dict] = {}
        with path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                species = row["species"]
                role = row["role"]
                roles[species] = {
                    "role": role,
                    "has_gff3": row["has_gff3"] == "yes",
                    "is_source": role in ("source", "both"),
                    "is_target": role in ("target", "both"),
                }
        if not roles:
            sys.exit(f"ERROR: no species rows in {path}")
        return roles

    @staticmethod
    def _read_top_sources(paths: Sequence[Path]) -> Dict[str, List[str]]:
        """{target: [source, …] in rank order} from <target>.top_sources.csv."""
        picks: Dict[str, List[str]] = {}
        for path in paths:
            target = path.name[: -len(".top_sources.csv")]
            with path.open(encoding="utf-8") as fh:
                picks[target] = [r["source"] for r in csv.DictReader(fh)]
        return picks

    # ── indexed views ──────────────────────────────────────────────────────────
    def annotation_index(self) -> Dict[tuple, dict]:
        """{(species, token, transcript_type): row} over the input annotations.

        Target-side rows carry a `_target` suffix on the token so a `role: both`
        species does not clobber its own source row (conf/modules.config adds it for
        exactly that reason); it is stripped here and the row's `role` column keeps
        the two apart.
        """
        index: Dict[tuple, dict] = {}
        for row in self.input_rows:
            parsed = split_sample(row["Sample"])
            if parsed is None:
                sys.exit(f"ERROR: unparseable Sample '{row['Sample']}' in {row['__file']}")
            species, token, ttype = parsed
            token = token[: -len("_target")] if token.endswith("_target") else token
            index[(species, row.get("role", "source"), token, ttype)] = row
        return index

    def gene_index(self) -> Dict[tuple, dict]:
        index: Dict[tuple, dict] = {}
        for row in self.gene_rows:
            parsed = split_sample(row["Sample"])
            if parsed is None:
                sys.exit(f"ERROR: unparseable Sample '{row['Sample']}' in {row['__file']}")
            species, token, category = parsed
            token = token[: -len("_target")] if token.endswith("_target") else token
            index[(species, row.get("role", "source"), token, category)] = row
        return index

    def projection_index(self) -> Dict[tuple, List[dict]]:
        """{(target, model, token): [rows]} — one row per transcript type, so the
        model count for a class is the sum over its rows."""
        index: Dict[tuple, List[dict]] = defaultdict(list)
        for row in self.projection_rows:
            parsed = split_sample(row["Sample"])
            if parsed is None:
                sys.exit(f"ERROR: unparseable Sample '{row['Sample']}' in {row['__file']}")
            prefix, token, _ttype = parsed
            if ".from_" not in prefix:
                sys.exit(f"ERROR: no '.from_' in projection Sample '{row['Sample']}'")
            target, model = prefix.split(".from_", 1)
            index[(target, model, token)].append(row)
        return index

    def td2_index(self) -> Dict[tuple, dict]:
        """{(sample, stage): row}; sample is the species for stage=input and
        `<target>.allModels` for stage=projected (one TD2 pass per target since the
        single-mapping redesign)."""
        return {(row["Sample"], row["stage"]): row for row in self.td2_rows}

    def alignment_index(self) -> Dict[tuple, dict]:
        """{(target, token): row} from `<target>.allModels.<token>`."""
        index: Dict[tuple, dict] = {}
        for row in self.align_rows:
            parts = row["Sample"].split(".")
            if len(parts) < 3:
                continue
            index[(".".join(parts[:-2]), parts[-1])] = row
        return index

    def accuracy(self) -> Dict[tuple, dict]:
        """{(target, source, feature_type, decoy): {sn, pr, f1}} at Transcript level."""
        acc: Dict[tuple, dict] = {}
        for path in self.stats_files:
            ident = parse_identity(str(path))
            if ident is None:
                print(f"WARNING: cannot parse identity from '{path}', skipping",
                      file=sys.stderr)
                continue
            sn_pr = transcript_sn_pr(str(path))
            if sn_pr is None:
                print(f"WARNING: no Transcript level in '{path}', skipping",
                      file=sys.stderr)
                continue
            sn, pr = sn_pr
            acc[(ident.target, ident.source, ident.feature_type, ident.decoy)] = {
                "sn": sn, "pr": pr, "f1": f1(sn, pr),
            }
        return acc


# ── section builders ────────────────────────────────────────────────────────────
def section_overview(inp: Inputs, args, summary: dict) -> None:
    roles = inp.roles
    sources = sorted(s for s, r in roles.items() if r["is_source"])
    targets = sorted(s for s, r in roles.items() if r["is_target"])
    annotated = [t for t in targets if roles[t]["has_gff3"]]
    self_pairs = [s for s in sources if s in targets]

    feature_types = ["lnc_RNA"] + (["mRNA"] if args.include_mrna else [])

    metrics = [
        ("Species in samplesheet", len(roles)),
        ("… role source (donor only)", sum(1 for r in roles.values() if r["role"] == "source")),
        ("… role target (annotated only)", sum(1 for r in roles.values() if r["role"] == "target")),
        ("… role both (donor and target)", sum(1 for r in roles.values() if r["role"] == "both")),
        ("Sources (S)", len(sources)),
        ("Targets (T)", len(targets)),
        ("Targets benchmarked (with GFF3)", len(annotated)),
        ("Source × target pairs", len(sources) * len(targets)),
        ("… of which self-pairs", len(self_pairs)),
        ("Feature types transferred", ", ".join(feature_types)),
        ("Decoy track (--include_decoy)", args.include_decoy),
        ("Top-N consensus (--top_consensus_n)", args.top_n),
        ("PSAURON cutoff (--td2_psauron_min)", args.psauron_min),
    ]
    write_tsv("run_overview_mqc.tsv", ["Metric", "Value"], metrics)

    summary["overview"] = {
        "n_species": len(roles),
        "n_source_only": sum(1 for r in roles.values() if r["role"] == "source"),
        "n_target_only": sum(1 for r in roles.values() if r["role"] == "target"),
        "n_both": sum(1 for r in roles.values() if r["role"] == "both"),
        "n_sources": len(sources),
        "n_targets": len(targets),
        "n_targets_annotated": len(annotated),
        "n_pairs": len(sources) * len(targets),
        "n_self_pairs": len(self_pairs),
        "feature_types": feature_types,
        "include_mrna": args.include_mrna,
        "include_decoy": args.include_decoy,
        "top_consensus_n": args.top_n,
        "td2_psauron_min": args.psauron_min,
    }
    summary["sources"] = sources
    summary["targets"] = targets


def _own_annotation(species: str, role_info: dict, ann: Dict[tuple, dict],
                    genes: Dict[tuple, dict]) -> dict:
    """A species' own input-annotation numbers, read from whichever side ran the
    stats: the source side if it donates, else the target side. A `role: both`
    species has both, and they are the same annotation."""
    side = "source" if role_info["is_source"] else "target"
    raw = ann.get((species, side, "1_raw", "lnc_RNA"))
    kept = ann.get((species, side, "2_lncRNA", "lnc_RNA"))
    gene_total = None
    for category in ("coding", "pseudogene", "non_coding"):
        row = genes.get((species, side, "1_raw", category))
        if row is not None:
            gene_total = add(gene_total or 0.0, num(row.get("n_genes")))
    return {
        "raw_row": raw,
        "kept_row": kept,
        "genes_total": gene_total,
    }


def section_species(inp: Inputs, args, summary: dict) -> None:
    ann, genes, td2 = inp.annotation_index(), inp.gene_index(), inp.td2_index()

    rank1 = defaultdict(int)
    in_top = defaultdict(int)
    for picks in inp.top_sources.values():
        for i, source in enumerate(picks):
            in_top[source] += 1
            if i == 0:
                rank1[source] += 1

    header = ["Sample", "role", "has_gff3", "is_source", "is_target", "benchmarked",
              "input_genes_total", "input_lncRNA_transcripts_raw", "spliced_kept",
              "td2_coding_dropped", "lncRNA_donated", "times_rank1", "times_in_top_n"]
    rows, records = [], {}
    for species in sorted(inp.roles):
        info = inp.roles[species]
        own = _own_annotation(species, info, ann, genes)
        td2_row = td2.get((species, "input"))
        record = {
            "role": info["role"],
            "has_gff3": info["has_gff3"],
            "is_source": info["is_source"],
            "is_target": info["is_target"],
            "benchmarked": info["is_target"] and info["has_gff3"],
            "input_genes_total": own["genes_total"],
            "input_lncRNA_transcripts_raw": num((own["raw_row"] or {}).get("n_transcripts")),
            "spliced_kept": num((td2_row or {}).get("n_in")),
            "td2_coding_dropped": num((td2_row or {}).get("n_coding")),
            "lncRNA_donated": num((td2_row or {}).get("n_kept")),
            "times_rank1": rank1.get(species, 0) if info["is_source"] else None,
            "times_in_top_n": in_top.get(species, 0) if info["is_source"] else None,
        }
        records[species] = record
        rows.append([species] + [record[k] for k in header[1:]])

    write_tsv("run_species_mqc.tsv", header, rows)
    summary["species"] = records


def section_source_funnel(inp: Inputs, args, summary: dict) -> None:
    ann, td2 = inp.annotation_index(), inp.td2_index()
    sources = [s for s, r in inp.roles.items() if r["is_source"]]

    header = ["Sample", "raw_lncRNA_transcripts", "dropped_single_exon", "dropped_coding",
              "kept", "pct_kept", "n_genes_raw", "n_genes_kept",
              "mean_transcript_length_kept", "mean_spliced_length_kept",
              "mean_exon_length_kept", "exons_per_transcript_kept"]
    rows, plot, records = [], {}, {}
    for species in sorted(sources):
        raw_row = ann.get((species, "source", "1_raw", "lnc_RNA"))
        kept_row = ann.get((species, "source", "2_lncRNA", "lnc_RNA"))
        td2_row = td2.get((species, "input"))

        raw = num((raw_row or {}).get("n_transcripts"))
        spliced = num((td2_row or {}).get("n_in"))
        coding = num((td2_row or {}).get("n_coding"))
        kept = num((td2_row or {}).get("n_kept"))
        single_exon = None if (raw is None or spliced is None) else raw - spliced

        # The GFF side must agree with the TD2 side: FILTER_LNC_GFF subsets the
        # filtered GFF3 by exactly the coding ids TD2 dropped. A mismatch means one
        # of the two stopped describing the same set — worth saying out loud.
        gff_kept = num((kept_row or {}).get("n_transcripts"))
        if kept is not None and gff_kept is not None and abs(kept - gff_kept) > 0.5:
            print(f"WARNING: {species}: TD2 kept {kept:.0f} lncRNA but the filtered "
                  f"GFF3 has {gff_kept:.0f} — the two filters have diverged",
                  file=sys.stderr)

        n_exons = num((kept_row or {}).get("n_exons"))
        record = {
            "raw_lncRNA_transcripts": raw,
            "dropped_single_exon": single_exon,
            "dropped_coding": coding,
            "kept": kept,
            "pct_kept": pct(kept, raw),
            "n_genes_raw": num((raw_row or {}).get("n_genes")),
            "n_genes_kept": num((kept_row or {}).get("n_genes")),
            "mean_transcript_length_kept": num((kept_row or {}).get("mean_transcript_length")),
            "mean_spliced_length_kept": num((kept_row or {}).get("mean_spliced_length")),
            "mean_exon_length_kept": num((kept_row or {}).get("mean_exon_length")),
            "exons_per_transcript_kept": (
                None if (n_exons is None or not gff_kept) else n_exons / gff_kept
            ),
        }
        records[species] = record
        rows.append([species] + [record[k] for k in header[1:]])

        # Stacked bars only make sense when the whole funnel resolved.
        if None not in (kept, single_exon, coding):
            plot[species] = {
                "kept": kept,
                "dropped_single_exon": single_exon,
                "dropped_coding": coding,
            }

    write_tsv("run_source_funnel_mqc.tsv", header, rows)
    write_json("run_source_funnel_mqc.json", {
        "id": "run_source_funnel_plot",
        "section_name": "Source lncRNA — what survived filtering",
        "description": (
            "Per source species: lncRNA transcripts donated (kept) and the two "
            "reasons the rest were dropped — single-exon (FOMO transfers spliced "
            "models only) and coding potential (a complete ORF with PSAURON score "
            "at or above the cutoff). Bars sum to the raw lncRNA count."
        ),
        "plot_type": "bargraph",
        "pconfig": {
            "id": "run_source_funnel_bargraph",
            "title": "FOMO: source lncRNA filtering funnel",
            "ylab": "transcripts",
            "cpswitch_counts_label": "Transcripts",
        },
        "data": plot,
    })
    summary["source_funnel"] = records


def section_targets(inp: Inputs, args, summary: dict) -> None:
    ann, projections, align = inp.annotation_index(), inp.projection_index(), inp.alignment_index()
    targets = sorted(s for s, r in inp.roles.items() if r["is_target"])
    # Only the ENABLED tracks can gain models. A type FOMO does not transfer gains a
    # known zero; an enabled type with no projection row is a genuine unknown (NA).
    enabled = {"lnc_RNA"} | ({"mRNA"} if args.include_mrna else set())

    def models(target: str, model: str, token: str) -> Optional[float]:
        rows = projections.get((target, model, token))
        if not rows:
            return None
        return add(*[num(r.get("n_transcripts")) or 0.0 for r in rows])

    header = ["Sample", "target", "transcript_type", "before",
              "added_allModels_collapsed", "added_top3_collapsed",
              "after_allModels", "after_top3",
              "pct_increase_allModels", "pct_increase_top3", "reads_mapped_percent"]
    rows, records = [], defaultdict(dict)
    bar_all, bar_top = {}, {}

    for target in targets:
        annotated = inp.roles[target]["has_gff3"]
        # Types present in the target's own annotation, plus the classes that
        # received projections (an un-annotated target only has the latter).
        types = sorted({
            key[3] for key in ann
            if key[0] == target and key[1] == "target" and key[2] == "1_raw"
        })
        for ttype in TYPE_TO_TOKEN:
            if ttype not in types and models(target, "allModels_collapsed", TYPE_TO_TOKEN[ttype]):
                types.append(ttype)
        types = sorted(set(types))

        for ttype in types:
            token = TYPE_TO_TOKEN.get(ttype) if ttype in enabled else None
            before = num((ann.get((target, "target", "1_raw", ttype)) or {}).get("n_transcripts")) \
                if annotated else None
            if token is None:
                added_all = added_top = 0.0
            else:
                added_all = models(target, "allModels_collapsed", token)
                added_top = models(target, "top3_collapsed", token)
            mapped = None
            if token:
                mapped = num((align.get((target, token)) or {}).get("reads_mapped_percent"))

            record = {
                "target": target,
                "transcript_type": ttype,
                "before": before,
                "added_allModels_collapsed": added_all,
                "added_top3_collapsed": added_top,
                "after_allModels": add(before, added_all),
                "after_top3": add(before, added_top),
                "pct_increase_allModels": pct(added_all, before),
                "pct_increase_top3": pct(added_top, before),
                "reads_mapped_percent": mapped,
            }
            records[target][ttype] = record
            rows.append([f"{target}.{ttype}"] + [record[k] for k in header[1:]])

        # The bargraph is the headline: existing lncRNA vs what FOMO added.
        lnc = records[target].get("lnc_RNA", {})
        existing = lnc.get("before") or 0.0
        if lnc.get("added_allModels_collapsed") is not None:
            bar_all[target] = {"existing_lncRNA": existing,
                               "new_candidates": lnc["added_allModels_collapsed"]}
        if lnc.get("added_top3_collapsed") is not None:
            bar_top[target] = {"existing_lncRNA": existing,
                               "new_candidates": lnc["added_top3_collapsed"]}

    write_tsv("run_target_types_mqc.tsv", header, rows)
    write_json("run_target_before_after_mqc.json", {
        "id": "run_target_before_after",
        "section_name": "Targets — lncRNA before and after FOMO",
        "description": (
            "Per target: lncRNA already annotated versus candidate lncRNA added by "
            "FOMO, for each collapsed aggregate (switch datasets above the plot). "
            "`allModels` pools every source; `top3` pools only the best-scoring ones, "
            "so it is a subset. A target with no annotation of its own starts at zero."
        ),
        "plot_type": "bargraph",
        "pconfig": {
            "id": "run_target_before_after_bargraph",
            "title": "FOMO: target lncRNA before and after",
            "ylab": "transcripts",
            "data_labels": [
                {"name": "allModels (every source)", "ylab": "transcripts"},
                {"name": "top3 (best sources)", "ylab": "transcripts"},
            ],
        },
        "data": [bar_all, bar_top],
    })
    summary["targets_before_after"] = {t: dict(v) for t, v in records.items()}


def section_heatmap(inp: Inputs, args, summary: dict) -> None:
    acc = inp.accuracy()
    sources = sorted(s for s, r in inp.roles.items() if r["is_source"])
    targets = sorted(s for s, r in inp.roles.items() if r["is_target"])

    # One section PER METRIC, not one plot with a dataset switcher: MultiQC 1.35
    # validates a custom-content heatmap's rows as scalars, so handing it a list of
    # matrices makes the whole report fail with "No analysis results found" (verified
    # against the 1.35 container). Bargraphs DO accept multi-dataset — heatmaps do not.
    METRICS = [
        ("f1", "F1", "F1 (the harmonic mean of the two below, and the score the "
                     "top-N ranking uses)"),
        ("sn", "Sensitivity", "Sensitivity — how much of the target's real annotation "
                              "the projection recovered"),
        ("pr", "Precision", "Precision — how much of the projection is real "
                            "annotation rather than noise"),
    ]
    shared = (
        "lncRNA transcript-level gffcompare accuracy for each source → target "
        "projection. A blank cell means that pair was not benchmarked — the target "
        "has no annotation of its own to compare against. The diagonal, where a "
        "species is projected onto itself, is the ceiling control: it is what a "
        "perfect donor would score, and it is excluded from every ranking."
    )

    records: Dict[str, dict] = {}
    for metric, label, blurb in METRICS:
        matrix = []
        for source in sources:
            row = []
            for target in targets:
                cell = acc.get((target, source, "lnc_RNA", False))
                value = None if cell is None else round(cell[metric], 2)
                row.append(value)
                if metric == "f1":
                    records.setdefault(source, {})[target] = value
            matrix.append(row)

        write_json(f"run_accuracy_heatmap_{metric}_mqc.json", {
            "id": f"run_accuracy_heatmap_{metric}",
            "section_name": f"Accuracy {label} — every source against every target",
            "description": f"{blurb}. {shared}",
            "plot_type": "heatmap",
            "xcats": targets,
            "ycats": sources,
            "pconfig": {
                "id": f"run_accuracy_heatmap_{metric}_plot",
                "title": f"FOMO: source × target lncRNA {label}",
                "xlab": "target",
                "ylab": "source",
                "min": 0,
                "max": 100,
            },
            "data": matrix,
        })
    summary["accuracy_f1"] = records


def section_aggregates(inp: Inputs, args, summary: dict) -> None:
    acc, projections = inp.accuracy(), inp.projection_index()
    targets = sorted(s for s, r in inp.roles.items() if r["is_target"])
    feature_types = ["lnc_RNA"] + (["mRNA"] if args.include_mrna else [])

    def models(target: str, model: str, ft: str, decoy: bool) -> Optional[float]:
        rows = projections.get((target, model, PROJECTION_TOKEN[(ft, decoy)]))
        if not rows:
            return None
        return add(*[num(r.get("n_transcripts")) or 0.0 for r in rows])

    header = ["Sample", "target", "feature_type", "model", "best_source", "sn", "pr", "f1",
              "n_models", "delta_f1_vs_best_single"]
    decoy_header = ["decoy_models", "decoy_f1", "est_fdr_pct"]
    any_decoy = args.include_decoy
    if any_decoy:
        header = header + decoy_header

    rows, plot, records = [], {}, defaultdict(dict)
    for target in targets:
        for ft in feature_types:
            # Baseline: the best real, non-self donor. Self is excluded for the same
            # reason consensus_top.nf excludes it from the ranking pool — a species
            # against itself is a near-perfect copy and would flatter nothing else.
            candidates = {
                source: cell["f1"]
                for (tgt, source, feat, decoy), cell in acc.items()
                if tgt == target and feat == ft and not decoy
                and source not in AGGREGATE_IDS and source != target
            }
            best_source = max(candidates, key=lambda s: (candidates[s], s), default=None)
            best_f1 = candidates.get(best_source) if best_source else None

            # (row label, the model id whose stats to read). The label is a category
            # name in the bargraph, so it must NOT carry the species — otherwise every
            # target contributes a differently-named bar and nothing is comparable.
            # Which species won goes in its own column instead.
            entries = []
            if best_source:
                entries.append(("best_single_source", best_source))
            entries += [(m, m) for m in AGGREGATE_IDS]

            for label, model in entries:
                cell = acc.get((target, model, ft, False))
                if cell is None:
                    continue
                real_models = models(target, model, ft, False)
                decoy_models = models(target, model, ft, True) if any_decoy else None
                decoy_cell = acc.get((target, model, ft, True)) if any_decoy else None
                record = {
                    "target": target,
                    "feature_type": ft,
                    "model": label,
                    "best_source": best_source,
                    "sn": cell["sn"],
                    "pr": cell["pr"],
                    "f1": cell["f1"],
                    "n_models": real_models,
                    "delta_f1_vs_best_single": (
                        None if best_f1 is None else cell["f1"] - best_f1
                    ),
                }
                if any_decoy:
                    record["decoy_models"] = decoy_models
                    record["decoy_f1"] = None if decoy_cell is None else decoy_cell["f1"]
                    record["est_fdr_pct"] = pct(decoy_models, real_models)
                records[target][f"{ft}.{label}"] = record
                rows.append([f"{target}.{ft}.{label}"] + [record[k] for k in header[1:]])

                if ft == "lnc_RNA":
                    plot.setdefault(target, {})[label] = round(cell["f1"], 2)

    write_tsv("run_aggregate_accuracy_mqc.tsv", header, rows)
    write_json("run_aggregate_accuracy_mqc.json", {
        "id": "run_aggregate_accuracy_plot",
        "section_name": "Aggregate models — lncRNA F1 per target",
        "description": (
            "lncRNA transcript-level F1 per target for the best single donor and for "
            "each aggregate model, side by side. This is the consensus question: does "
            "pooling every source (allModels), or only the top-N (top3), beat the best "
            "single donor — and does the gffcompare collapse help or hurt."
        ),
        "plot_type": "bargraph",
        "pconfig": {
            "id": "run_aggregate_accuracy_bargraph",
            "title": "FOMO: aggregate model accuracy (lncRNA, transcript level)",
            "ylab": "F1",
            "stacking": None,
            "cpswitch": False,
        },
        "data": plot,
    })
    summary["aggregates"] = {t: dict(v) for t, v in records.items()}


def section_donors(inp: Inputs, args, summary: dict) -> None:
    acc, projections = inp.accuracy(), inp.projection_index()
    sources = sorted(s for s, r in inp.roles.items() if r["is_source"])

    rank1, in_top = defaultdict(int), defaultdict(int)
    for picks in inp.top_sources.values():
        for i, source in enumerate(picks):
            in_top[source] += 1
            if i == 0:
                rank1[source] += 1

    header = ["Sample", "mean_f1", "median_f1", "best_f1", "n_targets_scored",
              "n_rank1", "n_in_top_n", "mean_models_projected"]
    rows, records = [], {}
    for source in sources:
        # Self-pairs excluded everywhere in this table (see consensus_top.nf).
        scores = [
            cell["f1"] for (target, src, ft, decoy), cell in acc.items()
            if src == source and ft == "lnc_RNA" and not decoy and target != source
        ]
        counts = [
            add(*[num(r.get("n_transcripts")) or 0.0 for r in rows_])
            for (target, model, token), rows_ in projections.items()
            if model == source and token == "1_lncRNA" and target != source
        ]
        counts = [c for c in counts if c is not None]
        record = {
            "mean_f1": statistics.fmean(scores) if scores else None,
            "median_f1": statistics.median(scores) if scores else None,
            "best_f1": max(scores) if scores else None,
            "n_targets_scored": len(scores),
            "n_rank1": rank1.get(source, 0),
            "n_in_top_n": in_top.get(source, 0),
            "mean_models_projected": statistics.fmean(counts) if counts else None,
        }
        records[source] = record
        rows.append([source] + [record[k] for k in header[1:]])

    write_tsv("run_donor_ranking_mqc.tsv", header, rows)
    write_json("run_donor_picks_mqc.json", {
        "id": "run_donor_picks",
        "section_name": "Donors — how often each source was picked",
        "description": (
            "How many targets picked each source into their top-N consensus, split by "
            "whether it ranked first. A source that is never picked contributes "
            "nothing to any `top3` model."
        ),
        "plot_type": "bargraph",
        "pconfig": {
            "id": "run_donor_picks_bargraph",
            "title": "FOMO: top-N donor selections per species",
            "ylab": "targets",
            "cpswitch_counts_label": "Targets",
        },
        "data": {
            s: {"picked_rank1": rank1.get(s, 0),
                "picked_rank2plus": in_top.get(s, 0) - rank1.get(s, 0)}
            for s in sources
        },
    })
    summary["donors"] = records

    # The actual selection per target, in rank order.
    width = max([len(p) for p in inp.top_sources.values()] or [0]) or args.top_n
    top_header = ["Sample"] + [f"rank_{i + 1}" for i in range(width)]
    top_rows = []
    for target in sorted(inp.top_sources):
        picks = inp.top_sources[target]
        top_rows.append([target] + [picks[i] if i < len(picks) else NA for i in range(width)])
    write_tsv("run_top_sources_mqc.tsv", top_header, top_rows)
    summary["top_sources"] = inp.top_sources


# ── entry point ─────────────────────────────────────────────────────────────────
def bool_arg(value: str) -> bool:
    return str(value).lower() in ("true", "1", "yes")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mqc-dir", required=True, type=Path,
                   help="directory holding the staged pipeline-wide MultiQC inputs")
    p.add_argument("--roles", required=True, type=Path,
                   help="species_roles.csv (species,role,has_gff3)")
    p.add_argument("--top-sources", nargs="*", type=Path, default=[],
                   help="<target>.top_sources.csv files (none if no target is annotated)")
    p.add_argument("--include-mrna", type=bool_arg, default=False)
    p.add_argument("--include-decoy", type=bool_arg, default=False)
    p.add_argument("--top-n", type=int, default=3)
    p.add_argument("--psauron-min", type=float, default=0.5)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not args.mqc_dir.is_dir():
        sys.exit(f"ERROR: --mqc-dir '{args.mqc_dir}' is not a directory")

    inp = Inputs(args.mqc_dir, args.roles, args.top_sources)
    summary: dict = {}

    section_overview(inp, args, summary)
    section_species(inp, args, summary)
    section_source_funnel(inp, args, summary)
    section_targets(inp, args, summary)
    section_heatmap(inp, args, summary)
    section_aggregates(inp, args, summary)
    section_donors(inp, args, summary)

    write_json("run_summary.json", summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
