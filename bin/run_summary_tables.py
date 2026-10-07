#!/usr/bin/env python3
"""
Digest a whole FOMO run into run-level MultiQC sections.

The per-target reports answer "how did this target do". This answers "what did this
run do": species per role, how much annotation went in and how much survived
filtering, how much each target gained, and which donors are worth using.

Everything is derived from artefacts the pipeline already produces — the
`*_mqc.tsv` adapter tables (including the curation tables, plans/21) and the
gffcompare `.stats` files — so this script adds no new measurement, only
cross-target arithmetic. Identity
comes from the `Sample` column (the ext.sample_name ladders in conf/modules.config
encode target/source/class) and from the stats filename grammar in fomo_stats.py.
Nothing here parses a GFF or a BAM.

Outputs, in report order (see assets/multiqc/run_summary_main.yml):
  run_overview_mqc.tsv            run overview: roles, pair counts, params
  run_species_mqc.tsv             per species: role + its own annotation input
  run_source_funnel_mqc.{tsv,json} per source: raw → stranded → spliced → longest → non-coding funnel
  run_target_types_mqc.tsv        per (target × transcript type): genes before / curated / after
  run_target_before_after_mqc.json  the lncRNA genes, as a stacked bargraph
  run_accuracy_heatmap_{f1,sn,pr}_mqc.json  source × target accuracy matrices
  run_aggregate_accuracy_mqc.{tsv,json}  allModels (every source, raw) vs best single donor
  run_donor_ranking_mqc.tsv       per source: F1 across targets
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

    def __init__(self, mqc_dir: Path, roles_csv: Path):
        files = sorted(p for p in mqc_dir.rglob("*") if p.is_file())
        self._assert_unique(files)

        def pick(suffix: str) -> List[Path]:
            return [p for p in files if p.name.endswith(suffix)]

        # role/kind/transcript-type tables of the INPUT annotations (source + target)
        self.input_rows = read_tsvs(pick("_gffstats_input_mqc.tsv"))
        self.gene_rows = read_tsvs(pick("_gffstats_genes_mqc.tsv"))
        self.projection_rows = read_tsvs(pick("_gffstats_projection_mqc.tsv"))
        self.td2_rows = read_tsvs(pick("_td2_coding_mqc.tsv"))
        # FILTER_TRANSCRIPT's own funnel (n_transcripts / n_spliced / n_selected),
        # one row per (source, feature_type), Sample = <species>.<feature_type>.
        self.filter_rows = read_tsvs(pick("_transcript_filter_mqc.tsv"))
        self.align_rows = read_tsvs(pick("_samtools_align_mqc.tsv"))
        self.stats_files = pick(".gffcompare.stats")
        # CURATE_MODELS (plans/21): one row per (target, track); Sample = <target>.<ft>[.decoy]
        self.curation_novel_rows = read_tsvs(pick("_curation_novel_mqc.tsv"))
        self.curation_observed_rows = read_tsvs(pick("_curation_observed_mqc.tsv"))

        self.roles = self._read_roles(roles_csv)

        print(
            "run_summary: staged "
            f"{len(self.input_rows)} input, {len(self.gene_rows)} gene, "
            f"{len(self.projection_rows)} projection, {len(self.td2_rows)} TD2, "
            f"{len(self.filter_rows)} transcript-filter, "
            f"{len(self.align_rows)} alignment rows; "
            f"{len(self.stats_files)} gffcompare stats; "
            f"{len(self.curation_novel_rows)} curation tracks",
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
    def _curation_key(sample: str) -> Optional[tuple]:
        """`<target>.<feature_type>[.decoy]` → (target, feature_type, decoy), suffixes
        stripped by length (species ids may contain anything but '.')."""
        decoy = sample.endswith(".decoy")
        rest = sample[: -len(".decoy")] if decoy else sample
        for ft in ("lnc_RNA", "mRNA"):
            if rest.endswith("." + ft):
                return rest[: -len(ft) - 1], ft, decoy
        return None

    def curation_index(self) -> Dict[tuple, dict]:
        """{(target, feature_type, decoy): {**novel row, **observed row}}."""
        index: Dict[tuple, dict] = {}
        for row in self.curation_observed_rows + self.curation_novel_rows:
            key = self._curation_key(row["Sample"])
            if key is None:
                sys.exit(f"ERROR: unparseable curation Sample '{row['Sample']}' in {row['__file']}")
            index.setdefault(key, {}).update(row)
        return index

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

    def filter_index(self) -> Dict[tuple, dict]:
        """{(species, feature_type): row}. Sample is `<species>.<feature_type>`; the
        species is recovered by stripping the KNOWN `.<feature>` suffix by length, not
        by splitting on '.', so a species id containing dots still round-trips."""
        index: Dict[tuple, dict] = {}
        for row in self.filter_rows:
            sample, feature = row["Sample"], row["feature"]
            suffix = f".{feature}"
            species = sample[:-len(suffix)] if sample.endswith(suffix) else sample
            index[(species, feature)] = row
        return index

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
    # Benchmarked = donor AND target with a gff3 (role 'both'): the gffcompare reference is the
    # species' filtered annotation, which only donors have. A pure target may carry a gff3 too
    # (plans/22) — it is then that target's curation reference, but not benchmarked.
    annotated = [t for t in targets if roles[t]["has_gff3"] and roles[t]["is_source"]]
    with_reference = [t for t in targets if roles[t]["has_gff3"]]
    # No self-projection (plans/23): a 'both' species is never projected onto itself.
    n_pairs = len(sources) * len(targets) - sum(1 for s in sources if s in targets)

    feature_types = ["lnc_RNA"] + (["mRNA"] if args.include_mrna else [])

    metrics = [
        ("Species in samplesheet", len(roles)),
        ("… role source (donor only)", sum(1 for r in roles.values() if r["role"] == "source")),
        ("… role target (projected onto, not a donor)", sum(1 for r in roles.values() if r["role"] == "target")),
        ("… role both (donor and target)", sum(1 for r in roles.values() if r["role"] == "both")),
        ("Sources (S)", len(sources)),
        ("Targets (T)", len(targets)),
        ("Targets with a reference annotation (curation)", len(with_reference)),
        ("Targets benchmarked (role both)", len(annotated)),
        ("Source × target pairs (no self-projection)", n_pairs),
        ("Feature types transferred", ", ".join(feature_types)),
        ("Decoy track (--include_decoy)", args.include_decoy),
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
        "n_targets_with_reference": len(with_reference),
        "n_pairs": n_pairs,
        "feature_types": feature_types,
        "include_mrna": args.include_mrna,
        "include_decoy": args.include_decoy,
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

    header = ["Sample", "role", "has_gff3", "is_source", "is_target", "benchmarked",
              "input_genes_total", "input_lncRNA_transcripts_raw", "spliced_kept",
              "td2_coding_dropped", "lncRNA_donated"]
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
            "benchmarked": info["is_target"] and info["is_source"] and info["has_gff3"],
            "input_genes_total": own["genes_total"],
            "input_lncRNA_transcripts_raw": num((own["raw_row"] or {}).get("n_transcripts")),
            "spliced_kept": num((td2_row or {}).get("n_in")),
            "td2_coding_dropped": num((td2_row or {}).get("n_coding")),
            "lncRNA_donated": num((td2_row or {}).get("n_kept")),
        }
        records[species] = record
        rows.append([species] + [record[k] for k in header[1:]])

    write_tsv("run_species_mqc.tsv", header, rows)
    summary["species"] = records


def section_source_funnel(inp: Inputs, args, summary: dict) -> None:
    ann, td2, filt = inp.annotation_index(), inp.td2_index(), inp.filter_index()
    sources = [s for s, r in inp.roles.items() if r["is_source"]]

    header = ["Sample", "raw_lncRNA_transcripts", "dropped_unstranded", "dropped_single_exon",
              "dropped_non_longest", "dropped_coding", "kept", "pct_kept", "n_genes_raw", "n_genes_kept",
              "mean_transcript_length_kept", "mean_spliced_length_kept",
              "mean_exon_length_kept", "exons_per_transcript_kept"]
    rows, plot, records = [], {}, {}
    for species in sorted(sources):
        raw_row = ann.get((species, "source", "1_raw", "lnc_RNA"))
        kept_row = ann.get((species, "source", "2_lncRNA", "lnc_RNA"))
        td2_row = td2.get((species, "input"))

        filt_row = filt.get((species, "lnc_RNA"))

        raw = num((raw_row or {}).get("n_transcripts"))
        coding = num((td2_row or {}).get("n_coding"))
        kept = num((td2_row or {}).get("n_kept"))
        # FILTER_TRANSCRIPT's funnel separates its three structural filters — strand
        # consistency (drops "?"/trans-spliced models), multi-exon, longest isoform
        # (plans/20). TD2's own n_in is what is LEFT after all three, so `raw - n_in`
        # alone would lump them together.
        n_tx = num((filt_row or {}).get("n_transcripts"))
        stranded = num((filt_row or {}).get("n_stranded"))
        spliced = num((filt_row or {}).get("n_spliced"))
        selected = num((filt_row or {}).get("n_selected"))
        unstranded = None if (n_tx is None or stranded is None) else n_tx - stranded
        single_exon = None if (stranded is None or spliced is None) else stranded - spliced
        non_longest = None if (spliced is None or selected is None) else spliced - selected

        # Each step must hand the next exactly what it says it kept. A mismatch means
        # two artefacts stopped describing the same set — worth saying out loud.
        td2_in = num((td2_row or {}).get("n_in"))
        if selected is not None and td2_in is not None and abs(selected - td2_in) > 0.5:
            print(f"WARNING: {species}: FILTER_TRANSCRIPT selected {selected:.0f} lncRNA "
                  f"but TD2 received {td2_in:.0f}", file=sys.stderr)
        if n_tx is not None and raw is not None and abs(n_tx - raw) > 0.5:
            print(f"WARNING: {species}: FILTER_TRANSCRIPT saw {n_tx:.0f} lncRNA but the "
                  f"raw GFF stats count {raw:.0f}", file=sys.stderr)
        # FILTER_ANNOTATION builds the filtered GFF3 from exactly the FASTA TD2 kept,
        # so this should be impossible now — kept as a regression check.
        gff_kept = num((kept_row or {}).get("n_transcripts"))
        if kept is not None and gff_kept is not None and abs(kept - gff_kept) > 0.5:
            print(f"WARNING: {species}: TD2 kept {kept:.0f} lncRNA but the filtered "
                  f"GFF3 has {gff_kept:.0f} — the two filters have diverged",
                  file=sys.stderr)

        n_exons = num((kept_row or {}).get("n_exons"))
        record = {
            "raw_lncRNA_transcripts": raw,
            "dropped_unstranded": unstranded,
            "dropped_single_exon": single_exon,
            "dropped_non_longest": non_longest,
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
        if None not in (kept, unstranded, single_exon, non_longest, coding):
            plot[species] = {
                "kept": kept,
                "dropped_unstranded": unstranded,
                "dropped_single_exon": single_exon,
                "dropped_non_longest": non_longest,
                "dropped_coding": coding,
            }

    write_tsv("run_source_funnel_mqc.tsv", header, rows)
    write_json("run_source_funnel_mqc.json", {
        "id": "run_source_funnel_plot",
        "section_name": "Source lncRNA — what survived filtering",
        "description": (
            "Per source species: lncRNA transcripts donated (kept) and the four "
            "reasons the rest were dropped — no consistent strand (\"?\" or "
            "trans-spliced across both strands; almost always zero), single-exon (FOMO transfers spliced "
            "models only), not the longest spliced isoform of its gene (one "
            "transcript per gene; zero with --longest_isoform false) and coding "
            "potential (a complete ORF with PSAURON score at or above the cutoff). "
            "Bars sum to the raw lncRNA count."
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
    """Per (target × transcript type): genes already annotated, curated genes FOMO adds
    (plans/21: intron chains shared by >= 2 species, none overlapping the reference), and
    the sum. Counted in GENES on both sides — a curated gene is one representative model,
    and the reference count is its gene count (isoforms would inflate 'before')."""
    ann, cur, align = inp.annotation_index(), inp.curation_index(), inp.alignment_index()
    targets = sorted(s for s, r in inp.roles.items() if r["is_target"])
    # Only the ENABLED tracks are curated. A type FOMO does not transfer gains a known
    # zero; an enabled type without a curation row is a genuine unknown (NA).
    enabled = {"lnc_RNA"} | ({"mRNA"} if args.include_mrna else set())

    def curated(target: str, ft: str, decoy: bool, key: str = "genes") -> Optional[float]:
        row = cur.get((target, ft, decoy))
        return None if row is None else num(row.get(key))

    header = ["Sample", "target", "transcript_type", "genes_before", "curated_genes",
              "genes_after", "pct_increase", "reads_mapped_percent"]
    if args.include_decoy:
        # Decoys are relocated per source into intergenic space, so a decoy curated gene is
        # a chance cross-species agreement. Normalised per multi-exon input model, because
        # decoy_cap makes the decoy and real inputs different sizes.
        header += ["curated_decoy_genes", "est_curated_fdr_pct"]
    rows, records = [], defaultdict(dict)
    bar = {}

    for target in targets:
        has_ref = inp.roles[target]["has_gff3"]
        types = {key[3] for key in ann if key[0] == target and key[1] == "target" and key[2] == "1_raw"}
        types |= {ft for (t, ft, d) in cur if t == target and not d}
        for ttype in sorted(types):
            before = num((ann.get((target, "target", "1_raw", ttype)) or {}).get("n_genes")) \
                if has_ref else None
            if ttype not in enabled:
                added = 0.0
            else:
                added = curated(target, ttype, False)
            token = TYPE_TO_TOKEN.get(ttype) if ttype in enabled else None
            mapped = num((align.get((target, token)) or {}).get("reads_mapped_percent")) if token else None
            record = {
                "target": target,
                "transcript_type": ttype,
                "genes_before": before,
                "curated_genes": added,
                "genes_after": add(before, added) if has_ref else added,
                "pct_increase": pct(added, before),
                "reads_mapped_percent": mapped,
            }
            if args.include_decoy:
                decoy = curated(target, ttype, True) if ttype in enabled else None
                real_in = curated(target, ttype, False, "multi_exon_models")
                decoy_in = curated(target, ttype, True, "multi_exon_models")
                rate_real = None if added is None or not real_in else added / real_in
                rate_decoy = None if decoy is None or not decoy_in else decoy / decoy_in
                record["curated_decoy_genes"] = decoy
                record["est_curated_fdr_pct"] = (
                    None if rate_real in (None, 0) or rate_decoy is None else 100.0 * rate_decoy / rate_real)
            records[target][ttype] = record
            rows.append([f"{target}.{ttype}"] + [record[k] for k in header[1:]])

        lnc = records[target].get("lnc_RNA", {})
        if lnc.get("curated_genes") is not None:
            bar[target] = {"existing_lncRNA_genes": lnc.get("genes_before") or 0.0,
                           "curated_lncRNA_genes": lnc["curated_genes"]}

    write_tsv("run_target_types_mqc.tsv", header, rows)
    write_json("run_target_before_after_mqc.json", {
        "id": "run_target_before_after",
        "section_name": "Targets — lncRNA genes before and after FOMO",
        "description": (
            "Per target: lncRNA genes already in its reference annotation versus the curated "
            "lncRNA genes FOMO adds — intron chains shared exactly by at least "
            "--curate_min_species source species, none overlapping a reference exon. A "
            "target without a reference annotation starts at zero."
        ),
        "plot_type": "bargraph",
        "pconfig": {
            "id": "run_target_before_after_bargraph",
            "title": "FOMO: target lncRNA genes before and after curation",
            "ylab": "genes",
        },
        "data": bar,
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
        ("f1", "F1", "F1 (the harmonic mean of the two below)"),
        ("sn", "Sensitivity", "Sensitivity — how much of the target's real annotation "
                              "the projection recovered"),
        ("pr", "Precision", "Precision — how much of the projection is real "
                            "annotation rather than noise"),
    ]
    shared = (
        "lncRNA transcript-level gffcompare accuracy for each source → target "
        "projection. A blank cell means that pair was not benchmarked — the target "
        "is not benchmarked (only role `both` targets are). The diagonal is always "
        "blank: a species is never projected onto itself (plans/23)."
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
            # Baseline: the best real donor. (`source != target` only matters for a pre-23
            # run that still carried self-projections — a near-perfect copy of itself.)
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
        "section_name": "allModels vs best single donor — lncRNA F1 per target",
        "description": (
            "lncRNA transcript-level F1 per target for the best single donor and for "
            "allModels (every source pooled, raw), side by side: does pooling every "
            "donor beat the best single one?"
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
    """Per source: lncRNA transcript F1 across the benchmarked targets, and how many
    models it projects per target (self-pairs never exist since plans/23)."""
    acc, projections = inp.accuracy(), inp.projection_index()
    sources = sorted(s for s, r in inp.roles.items() if r["is_source"])

    header = ["Sample", "mean_f1", "median_f1", "best_f1", "n_targets_scored",
              "mean_models_projected"]
    rows, records = [], {}
    for source in sources:
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
            "mean_models_projected": statistics.fmean(counts) if counts else None,
        }
        records[source] = record
        rows.append([source] + [record[k] for k in header[1:]])

    write_tsv("run_donor_ranking_mqc.tsv", header, rows)
    summary["donors"] = records


# ── entry point ─────────────────────────────────────────────────────────────────
def bool_arg(value: str) -> bool:
    return str(value).lower() in ("true", "1", "yes")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mqc-dir", required=True, type=Path,
                   help="directory holding the staged pipeline-wide MultiQC inputs")
    p.add_argument("--roles", required=True, type=Path,
                   help="species_roles.csv (species,role,has_gff3)")
    p.add_argument("--include-mrna", type=bool_arg, default=False)
    p.add_argument("--include-decoy", type=bool_arg, default=False)
    p.add_argument("--psauron-min", type=float, default=0.5)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not args.mqc_dir.is_dir():
        sys.exit(f"ERROR: --mqc-dir '{args.mqc_dir}' is not a directory")

    inp = Inputs(args.mqc_dir, args.roles)
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
