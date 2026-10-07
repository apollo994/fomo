# 16 — Run-level "general" MultiQC report

## Context

Since plans/14 an all-vs-all experiment is **one** `nextflow run` (`role: both`), and
since plans/15 there is one merged alignment per target. Reporting never caught up:
`REPORTING` emits **T** MultiQC reports under `targets/<target>/multiqc/`, source-side
sections broadcast into each. With ~20 species that is 20 near-identical HTML files and
**no artefact that answers "what did this run do, and did it work?"** — the run-level
numbers (how many species in which role, how much annotation went in, how much came out
per target, which donors are any good) exist only implicitly, spread across 20 reports
and ~1000 stat files.

This plan adds **one run-level report** — `results/summary/multiqc/fomo_run_summary.html`
— built from purpose-made cross-target aggregates, plus the machine-readable TSV/JSON
those aggregates are made of (useful directly for figures and downstream analysis). The
per-target path is left byte-identical.

**Explicitly not** a second dump of the raw union: at S=T=20 the union is ~1000 files /
~900 table rows. Feeding it to a second MultiQC would produce a *bigger* report than the
per-target ones, and the per-target scatter JSONs would collide on their constant section
id `gffcompare_accuracy_custom`. The general report is aggregate-only.

## Decisions

| Question | Decision |
|---|---|
| Where it lives | New top-level subworkflow `RUN_SUMMARY` (`subworkflows/local/run_summary.nf`), called from `workflows/fomo.nf`. **Not** nested inside `REPORTING` and not nested under anything — every `conf/modules.config` selector is a fully-qualified path (rule 1 in CLAUDE.md). |
| MultiQC instance | `include { MULTIQC as MULTIQC_RUN_SUMMARY }`. An alias is required: the publishDir selector `.*:MULTIQC$` is `$`-anchored so the alias escapes it, and `withName: 'FOMO:REPORTING:MULTIQC'` (unanchored) would otherwise regex-*find* inside a same-subworkflow alias. Separate subworkflow + alias makes both misses structural. |
| Where the numbers come from | The existing `*_mqc.tsv` adapters' output and the gffcompare `.stats` files — i.e. `PREPROCESSING/PROJECTION/BENCHMARKING/CONSENSUS_TOP.out.mqc_files`, which is **already** the complete pipeline-wide union. No new stat producer, no re-running any tool. |
| How identity is recovered | From the `Sample` column of the TSVs (the `ext.sample_name` ladders already encode target/source/class) and from the `.stats` filename grammar already used by `bin/select_top_sources.py` / `bin/gffcompare_accuracy_mqc.py`. No new naming contract. |
| One aggregator task | `RUN_SUMMARY_TABLES` — one task per run over `.collect()`, container `python:3.11` (stdlib only, same choice as `GFF_STATS_TO_MQC`; avoids a fourth place pinning the MultiQC container). |
| Output shape | TSVs for tables (deterministic column order, the established pattern), self-describing `*_mqc.json` custom content for plots (bargraph/heatmap), copying `bin/gffcompare_accuracy_mqc.py`'s JSON idiom. |
| Publishing | `results/summary/multiqc/fomo_run_summary.html` + `results/summary/tables/*.{tsv,json}`. Static paths — no `meta.target_id` segment (a target-scoped closure on a run-level process publishes to `targets/null/`). |
| "After" for the per-target before/after | **Both aggregates side by side**: `after_allModels = before + allModels.collapsed`, `after_top3 = before + top3.collapsed`, as separate columns. |
| Ranking/parse duplication | The identity regex + `AGGREGATE_IDS` + the stats parser exist in **two** bin scripts already; a third copy is not acceptable. Extract them into `bin/fomo_stats.py` and import from all three. |
| Deferred | Run-integrity (expected-vs-observed task counts) and software versions via `Channel.topic('versions')` — the topic is currently collected nowhere, so a versions section is free later. Not in this plan. |

## Files

**New**
- `bin/fomo_stats.py` — shared helpers: `FNAME_RE`, `AGGREGATE_IDS`, `parse_identity()` (now returning **target** too), `parse_accuracy()`, `f1()`.
- `bin/run_summary_tables.py` — the aggregator.
- `modules/local/run_summary_tables.nf` — `RUN_SUMMARY_TABLES`.
- `subworkflows/local/run_summary.nf` — `RUN_SUMMARY`.
- `assets/multiqc/run_summary_main.yml`, `assets/multiqc/run_summary_sections.yml`.

**Modified**
- `workflows/fomo.nf` — hoist the union into a named var; build the roles CSV; call `RUN_SUMMARY`.
- `subworkflows/local/consensus_top.nf` — add `top_sources = SELECT_TOP_SOURCES.out.csv` to `emit:`.
- `bin/gffcompare_accuracy_mqc.py`, `bin/select_top_sources.py` — import from `fomo_stats`.
- `conf/modules.config` — two new `withName` blocks + two new publishDir rules; `$`-anchor the existing `FOMO:REPORTING:MULTIQC` selector.
- `nextflow.config` + `nextflow_schema.json` — two new params (kept in sync or the run aborts: `failUnrecognisedParams = true`).
- `CLAUDE.md` — see the last section.

## 1. Plumbing (`workflows/fomo.nf`)

Pure refactor — DSL2 forks a channel read twice, so `REPORTING` is unaffected:

```groovy
    ch_all_mqc = PREPROCESSING.out.mqc_files
        .mix(PROJECTION.out.mqc_files)
        .mix(BENCHMARKING.out.mqc_files)
        .mix(CONSENSUS_TOP.out.mqc_files)

    REPORTING(ch_all_mqc, BENCHMARKING.out.stats.mix(CONSENSUS_TOP.out.stats), ch_target_ids)

    // Species roles are the one fact the stat files do not carry. Emit them as a
    // 3-column CSV from the samplesheet rows already read above; `seed` is written
    // first, `sort` makes the body deterministic.
    ch_roles = Channel.fromList(rows)
        .map { meta, _fa, gff3 -> "${meta.id},${meta.role},${gff3 ? 'yes' : 'no'}" }
        .collectFile(name: 'species_roles.csv', newLine: true, sort: true,
                     seed: 'species,role,has_gff3')

    RUN_SUMMARY(ch_all_mqc, CONSENSUS_TOP.out.top_sources, ch_roles)
```

`ch_target_ids` is today an inline expression in the `REPORTING(...)` call; hoist it too.

## 2. `subworkflows/local/run_summary.nf`

```groovy
workflow RUN_SUMMARY {
    take:
    ch_all_mqc      // Channel<tuple(meta, path|List<path>)> — the pipeline-wide union
    ch_top_sources  // Channel<tuple(meta, *.top_sources.csv)> × T_g
    ch_roles        // Channel<path> — species_roles.csv (single file)

    main:
    // flatMap, NOT map: GFF_STATS_TO_MQC items are Lists of two paths (the
    // transcript table and the gene-category table) — the same trap REPORTING
    // documents. Everything is then collected into ONE task; the script selects
    // files by suffix from its own work dir, so no per-suffix channel splitting.
    ch_files = ch_all_mqc.flatMap { _meta, p -> p instanceof List ? p : [p] }.collect()
    ch_top   = ch_top_sources.map { _meta, csv -> csv }.collect(sort: true)

    RUN_SUMMARY_TABLES(ch_files, ch_top, ch_roles)

    MULTIQC_RUN_SUMMARY(
        RUN_SUMMARY_TABLES.out.mqc.map { files -> tuple(
            [id: 'run_summary'],                       // NO target_id — deliberate
            files,
            [ file(params.multiqc_run_summary_main_config,     checkIfExists: true),
              file(params.multiqc_run_summary_sections_config, checkIfExists: true) ],
            [], [], []
        ) }
    )

    emit:
    report = MULTIQC_RUN_SUMMARY.out.report
    tables = RUN_SUMMARY_TABLES.out.tables
}
```

`ch_top_sources` may be **empty** (no annotated target ⇒ no ranking): `.collect()` on an
empty channel emits nothing and would deadlock the task. Use
`.collect().ifEmpty([])` for `ch_top` (and keep `ch_files` un-guarded — it is never
empty). `RUN_SUMMARY_TABLES` declares those inputs as `path(..., stageAs: '?/*')`-style
optional lists.

## 3. `RUN_SUMMARY_TABLES` (module + script)

Process: `label 'process_single'`, container `python:3.11`, inputs
`(path mqc_files, path top_sources, path roles_csv)`, outputs
`path("*_mqc.{tsv,json}"), emit: mqc` and `path("run_summary.json"), emit: tables`
(+ the usual `versions` topic line, `python3 --version`).

`ext.args` carries the params the numbers must be read against (the script never touches
`params`):

```groovy
withName: 'FOMO:RUN_SUMMARY:RUN_SUMMARY_TABLES$' {
    ext.args = { "--include-mrna ${params.include_mrna} --include-decoy ${params.include_decoy} " +
                 "--top-n ${params.top_consensus_n} --psauron-min ${params.td2_psauron_min}" }
}
```

### Inputs the script reads (all by suffix glob in the work dir)

| Glob | Gives |
|---|---|
| `*_gffstats_input_mqc.tsv` | per-species × kind × transcript_type counts and length stats. `role` column = `source`/`target`; target rows' `Sample` ends `_target` (`<t>.1_raw_target`, `<t>.2_lncRNA_target`) |
| `*_gffstats_genes_mqc.tsv` | gene-category counts (`coding`/`pseudogene`/`non_coding`) for the same inputs |
| `*_gffstats_projection_mqc.tsv` | projected model counts per `<target>.from_<source\|aggregate>.<class>` |
| `*_td2_coding_mqc.tsv` | `stage` ∈ {input, projected}: `n_in`, `n_coding`, `n_kept`, `pct_coding` |
| `*_samtools_align_mqc.tsv` | `reads_mapped_percent` per (target, class) |
| `*.gffcompare.stats` | Sn/Pr per level, identity from the filename `<target>.from_<source>.<ft>[.decoy].gffcompare.stats` |
| `*.top_sources.csv` | the actual top-N pick per target (header `source`, ≤N rows) |
| `species_roles.csv` | species, role, has_gff3 |

Two robustness rules, both worth a loud failure rather than a silent wrong number:
- **assert unique basenames** across the staged set (they are unique today — every
  target-scoped `ext.prefix` starts with `${meta.target_id}.` and source-side ones are
  species-keyed); duplicates mean a naming contract broke.
- an F1/count that cannot be resolved is written as `NA`, never `0`; a *row* that cannot
  be identified is a hard error naming the file.

### Outputs (sections, in report order)

1. **`run_overview_mqc.tsv`** — metric/value table (**how many species in which role**):
   `n_species`, `n_source_only`, `n_target_only`, `n_both`, `n_sources` (S), `n_targets`
   (T), `n_targets_annotated` (T_g), `n_pairs` (S·T), `n_self_pairs`, `feature_types`,
   `include_mrna`, `include_decoy`, `top_consensus_n`, `td2_psauron_min`.
2. **`run_species_mqc.tsv`** — one row per species: `role`, `has_gff3`, `is_source`,
   `is_target`, `benchmarked`, `input_genes_total`, `input_lncRNA_transcripts_raw`,
   `spliced_kept`, `td2_coding_dropped`, `lncRNA_donated`, `times_in_top_n`,
   `times_rank1`.
3. **`run_source_funnel_mqc.tsv` + `run_source_funnel_mqc.json`** (**genes extracted and
   filtered out, with summary stats**) — the extraction/filtering funnel per source with
   the drop attributed: `raw` → `dropped_single_exon` = `raw − td2.n_in` →
   `dropped_coding` = `td2.n_coding` → `kept` = `td2.n_kept`, plus `pct_kept`,
   `mean_transcript_length`, `mean_spliced_length`, `mean_exon_length`,
   `n_exons_per_transcript`, `n_genes` (raw vs kept). JSON is the stacked bargraph of
   `[kept, dropped_single_exon, dropped_coding]`.
   *Consistency check to log:* `kept` must equal the `kind=filtered` lncRNA
   `n_transcripts` from `gffstats_input` — `FILTER_LNC_GFF` applies exactly that drop.
4. **`run_target_types_mqc.tsv`** (**before/after gene-type count per target**, primary)
   — one row per (target × transcript_type): `before` (target reference, `kind=raw`),
   `added_allModels_collapsed`, `added_top3_collapsed`, `after_allModels`, `after_top3`,
   `pct_increase_allModels`, `pct_increase_top3`. Non-lncRNA types get `added = 0`.
   Un-annotated targets have `before = NA` and only the `added_*` columns — the
   structural consequence of no `gff3`, shown rather than hidden.
5. **`run_target_before_after_mqc.json`** (before/after at a glance) — bargraph, two
   switchable datasets (`allModels` / `top3` via `pconfig.data_labels`), each stacked
   `[existing lncRNA, new candidates]` per target. `reads_mapped_percent` rides along as
   an extra column in the table above, so a target that simply aligned badly is one
   glance away.
6. **`run_accuracy_heatmap_mqc.json`** — source × target heatmap, three switchable
   datasets **F1 / Sn / Pr** at `Transcript` level, `lnc_RNA`, non-decoy. Rows = source
   species (aggregates excluded — they are section 7), columns = target, diagonal =
   self-projection ceiling. Note the format check below.
7. **`run_aggregate_accuracy_mqc.tsv` + `.json`** — consensus gain: one row per
   (target × model ∈ `best_single_source`, `allModels_raw`, `allModels_collapsed`,
   `top3_raw`, `top3_collapsed`) with `Sn`, `Pr`, `F1`, `n_models`, and
   `delta_F1_vs_best_single`. JSON = grouped F1 bargraph per target. With
   `--include_decoy`, add `decoy_models`, `decoy_F1`, `est_fdr =
   decoy_models / real_models` columns (absent otherwise — sparse-safe, no rows, no
   section noise).
8. **`run_donor_ranking_mqc.tsv` + `run_top_sources_mqc.tsv`** — per source:
   `mean_f1`, `median_f1`, `n_targets_scored`, `n_rank1`, `n_in_top_n`,
   `mean_models_projected` (self-pairs excluded from every one of those, matching
   `consensus_top.nf`'s `meta.id != meta.target_id`); and per target the actual picks as
   `rank_1 … rank_N`. Plus a bargraph of times-picked per species.
9. **`run_summary.json`** — every number above in one machine-readable file (not a
   MultiQC input; published to `summary/tables/`).

**Format check before writing the heatmap:** MultiQC 1.35 custom-content heatmaps accept
either `data` as a list-of-lists with `xcats`/`ycats`, or a dict-of-dicts. Write the
dict-of-dicts form first and verify it renders in the test run; fall back to
matrix + `xcats`/`ycats` if not. Same "verify in the test run" applies to the
multi-dataset bargraphs — `pconfig.data_labels` is proven for the scatter
(`bin/gffcompare_accuracy_mqc.py`) but not yet for a bargraph here.

## 4. MultiQC config (`assets/multiqc/run_summary_*.yml`)

`run_summary_main.yml`: `report_title: "FOMO run summary"`, a `report_comment` stating
that self-projections are included in `allModels` and excluded from the top-N pool,
`skip_generalstats: true`, `report_section_order` numbering sections 1→9 top-down, and
`extra_fn_clean_exts` for the new `_mqc`/`.tsv` suffixes.

`run_summary_sections.yml`: one `custom_data` block + one `sp` pattern per **TSV**
(`*_run_overview_mqc.tsv`, `*_run_species_mqc.tsv`, …). The JSON plots are
self-describing and get **no** `sp` entry and **no** `custom_data` block — an id embedded
in the file plus an `sp` match renders the section twice (CLAUDE.md, learned the hard
way). Every new suffix is globally unique, so the per-target config can never route them
even though the two reports never share a task.

Params (both `hidden: true`, `format: file-path`, `exists: true`, defaults under
`${projectDir}/assets/multiqc/`): `multiqc_run_summary_main_config`,
`multiqc_run_summary_sections_config` — added to `nextflow.config` **and**
`nextflow_schema.json` together.

## 5. `conf/modules.config`

```groovy
withName: 'FOMO:RUN_SUMMARY:MULTIQC_RUN_SUMMARY$' {
    ext.prefix = 'fomo_run_summary'          // → fomo_run_summary.html
}
```

and, in the "Result publishing" section at the bottom (static paths, same `saveAs` that
drops `versions.yml`):

```groovy
withName: '.*:RUN_SUMMARY_TABLES$'   { publishDir = [ path: "${params.outdir}/summary/tables",  … ] }
withName: '.*:MULTIQC_RUN_SUMMARY$'  { publishDir = [ path: "${params.outdir}/summary/multiqc", … ] }
```

Defensive one-char change: anchor the existing `withName: 'FOMO:REPORTING:MULTIQC'` as
`'FOMO:REPORTING:MULTIQC$'` so a future same-subworkflow alias cannot inherit
`ext.prefix = 'multiqc_report'` by substring match.

## 6. Implementation order

1. `bin/fomo_stats.py`; repoint `gffcompare_accuracy_mqc.py` and `select_top_sources.py`
   at it (`sys.path.insert(0, str(Path(__file__).resolve().parent))` — `bin/` is staged
   as a directory on every executor, so sibling import is safe); `parse_identity` gains
   the target field. Re-run `-profile test` — the per-target reports must be unchanged.
2. `CONSENSUS_TOP` gains the `top_sources` emit; `fomo.nf` gains `ch_all_mqc`,
   `ch_roles`, `RUN_SUMMARY`; module + subworkflow + a *stub* script that only writes
   `run_overview_mqc.tsv`. Confirm the task runs once, the report lands in
   `results/summary/multiqc/`, and nothing appears under `targets/null/`.
3. Fill in the script section by section, in the numbered order above; re-run with
   `-resume` after each (only `RUN_SUMMARY_TABLES` + `MULTIQC_RUN_SUMMARY` re-execute).
4. `CLAUDE.md`.

## 7. Verification

```sh
nextflow run . -profile test,singularity --outdir results_summary        # T=3, S=4, F=1, D=1
nextflow run . -profile test,singularity --outdir results_summary_full \
    --include_mrna --include_decoy -resume                              # F=2, D=2 branches
```

`-profile test` is the right harness: 2 `both` rows (⇒ self-pairs on the diagonal), 1
pure `target` with an **empty** gff3 (⇒ no gffcompare, no top-N — the `NA` paths in
sections 4/6/7/8), and 2 pure sources.

Check, in order:
- `results_summary/summary/multiqc/fomo_run_summary.html` exists; `results_summary/targets/*/multiqc/multiqc_report.html` are **unchanged** vs. a pre-change run (diff the `multiqc_data/multiqc_data.json`, ignoring timestamps); nothing under `results_summary/targets/null/`.
- `run_overview_mqc.tsv`: `n_species 5`, `n_both 2`, `n_source_only 2`, `n_target_only 1`, `n_sources 4`, `n_targets 3`, `n_targets_annotated 2`, `n_pairs 12`, `n_self_pairs 2`.
- `run_species_mqc.tsv`: the un-annotated target has `has_gff3 no`, `is_source no`, `benchmarked no`, and `NA` (not 0) for every input-annotation column.
- Funnel arithmetic per source: `kept + dropped_single_exon + dropped_coding == raw`, and `kept` equals the `kind=filtered` lncRNA `n_transcripts` in `*_gffstats_input_mqc.tsv`.
- Before/after cross-check on one annotated target: `before` for `lnc_RNA` equals that target's `<t>.1_raw_target` row, `added_allModels_collapsed` equals its `from_allModels_collapsed.1_lncRNA` projection row, and `after = before + added`.
- Heatmap: 4 source rows × 3 target columns, the un-annotated target's column all `NA`, the two diagonal (self) cells the highest in their column — the ceiling control behaving as documented.
- Aggregate table: `top3_*` rows absent for the un-annotated target; `delta_F1_vs_best_single` finite for the other two.
- Donor table: no self-pair contributes to `mean_f1` (a `both` species' own column is skipped), and `n_in_top_n` totals match the `top_sources.csv` files.
- Second run: `mRNA` and `decoy_*` rows/columns appear; `est_fdr` column present; per-target reports still fine.
- `nextflow run . --help` renders the two new params (proves the schema is in sync).

## 8. `CLAUDE.md` updates

- "Reporting" stage in the architecture list: MultiQC report **per target plus one run-level summary**.
- Subworkflow map: new `RUN_SUMMARY` row (aggregates the union into run-level tables; one task per run).
- "Reporting conventions": `REPORTING` is no longer "the only place that calls MULTIQC" — `RUN_SUMMARY` calls the alias `MULTIQC_RUN_SUMMARY` (`RUN_SUMMARY_TABLES` runs on `python:3.11`, so the MultiQC container is still pinned in one place).
- Result-publishing tree: add the `summary/` sibling of `targets/`.
- Replace the "three places that must agree" note about `AGGREGATE_IDS` with **two**: `consensus_top.nf`'s channel filter and `bin/fomo_stats.py`.
- Note that the run-level report is aggregate-only by design, and why (union size + the constant `gffcompare_accuracy_custom` section id).
- Params & validation: the two new config params.

---

## As built — deviations from the plan above

Implemented and verified 2026-08-13 on `-profile test,crg`. Six differences, all
deliberate:

1. **The accuracy heatmap is three sections, not one with a dataset switcher.** The
   plan's format check found the answer: MultiQC 1.35 validates a custom-content
   heatmap's `rows` as scalars, so a list of matrices fails validation and the *whole*
   report comes out as "No analysis results found". Multi-dataset **bargraphs** do work
   (`pconfig.data_labels`), and `run_target_before_after` uses one. So the matrix ships
   as `run_accuracy_heatmap_{f1,sn,pr}_mqc.json`, three single-matrix sections.
2. **`est_fdr` → `est_fdr_pct`.** It is computed as a percentage, and on small data it
   legitimately exceeds 100% (relocated decoys don't overlap, so the collapse merges
   none of them, and the decoy track can be as large as the real one). The name now
   says so.
3. **`run_overview_mqc.tsv` uses prose metric labels** ("Species in samplesheet", not
   `n_species`) — it is the first thing a reader sees. The machine-readable keys
   (`n_species`, `n_targets_annotated`, …) are in `run_summary.json`, which is what a
   script should read.
4. **Three added columns**, all supersets of the plan: `feature_type` and `best_source`
   on the aggregate table (so mRNA rows appear under `--include_mrna`, and which donor
   won is explicit), and `best_f1` on the donor table. The aggregate **bargraph**
   categories use the bare `best_single_source` label — with the species name in the
   category, every target contributed a differently-named bar and nothing was
   comparable.
5. **No stub-script step.** Instead the aggregator was verified offline against a
   synthetic fixture covering both flag combinations, and its output rendered through
   the MultiQC 1.35 container, *before* wiring — cheaper and stricter than a stub, and
   it is what caught the heatmap format and the bargraph-category problems.
6. **The verification numbers in §7 were wrong.** They were taken from the header
   comment in `conf/test.config`, which claims "2 rows 'both' … T = 3, S = 4";
   `assets/samplesheet.csv` actually holds **3** species — 1 `both`, 1 pure `target`
   with no gff3, 1 pure `source` — so S=2, T=2, T_g=1, 4 pairs, 1 self-pair. Verified
   against those. The profile still exercises everything the plan wanted: a self-pair
   (the diagonal reads 92.92 vs the best real donor's 27.49) and an un-annotated target
   (blank heatmap column, `before = NA`, no `top3` rows).

Two adjacent defects found and **not** fixed (out of scope, worth a follow-up):
`conf/test.config`'s header comment is stale as above, and `nextflow.config` registers a
`test_single` profile whose `conf/test_single.config` does not exist — `-profile
test_single` therefore aborts.
