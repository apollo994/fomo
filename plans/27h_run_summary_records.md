# 27h — The run summary reads data records

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decision: U7. Review items:
P2, P4, D2, R5, and F5/F6/F8–F13/F16 of the Python review.

## Status: done (2026-10-15)

## Goal

`RUN_SUMMARY_TABLES` builds every table and `run_summary.json` from **records only** (27g) plus the
roles CSV. It no longer parses Sample labels (`1_raw`, `2_lncRNA`, `PROJECTION_TOKEN`),
filename suffixes or `.from_` grammars. The code is split into cached record indexes, pure
`compute_<section>` functions and `emit_<section>` writers. **Outputs stay byte-identical (U1).**

## Steps

- [x] **A1. Extend the regression fixture from 27a** (`tests/run_summary_tables/fixtures/{V1,V2}/`)
  with the records of the same runs (captured in 27g E1). (updated by D27g-5: every produced
  `*.data.json` of the V1/V2 gate runs is in `fomo_refactor/gate/27g/V{1,2}/work`, skip staged
  symlinks; the curation record carries `has_reference`; source-side raw gffstats records have
  no `track`.) The existing test (old inputs →
  golden outputs) keeps passing until A5.
- [x] **A2. Write the new input layer, `RecordIndex`.**
  - It loads every `*.data.json` once and builds the indexes: annotation stats by
    (species, side, kind, track); transcript funnel; TD2; alignment; accuracy by
    (target, source, track); curation by (target, track).
  - Use `functools.cached_property`, so no index is rebuilt; today `accuracy()` is rebuilt 3×
    and `annotation_index()` 3×.
  - **Target-side raw stats:** use the `side=target` record if present, otherwise the
    `side=source` record of the same species (identical file). This prepares 27i's removal of
    duplicate target stats (U5) without changing values.
- [x] **A3. Split every section into `compute_<section>(ctx) -> records` (pure) and
  `emit_<section>(records)`.**
  - Sections: overview, species, source_funnel, target_types, heatmap, aggregates, donors.
  - The funnel consistency checks become `check_funnel_consistency()`.
  - A `RunContext` holds sources, targets and enabled feature types, derived once.
  - Remove the 4 unused `args` parameters and the out-parameter `summary` dict mutation; `main`
    assembles the JSON.
- [x] **A4. Errors and naming.**
  - Raise `FomoError` (from `fomo_io`), not `sys.exit` inside the class; one policy for missing
    inputs.
  - Use the shared `fmt`/`pct` helpers and the shared NA policy.
  - Delete `PROJECTION_TOKEN`, `split_sample`, `_curation_key` and every hard-coded `"1_raw"`,
    `"2_lncRNA"` or `"1_lncRNA"`.
- [x] **A5. Switch `RUN_SUMMARY_TABLES`'s input.**
  - `run_summary.nf` passes the records (the union of the `data` emits) instead of `ch_all_mqc`.
  - Delete the old-input regression test, and rename the record-based one to be the regression
    test. Its expected outputs are the same golden `summary/tables`.
- [x] **A6. Remove the per-file ignores of `run_summary_tables.py`.** C901 ≤ 15 holds.
- [x] **A7. MultiQC_RUN_SUMMARY still renders.** Its inputs (`run_*_mqc.*`) are unchanged in name
  and content.

## Acceptance

- G1 (the fixture regression is byte-identical on V1 and V2), G2 (no ignores for this file), G3,
  and G4 identical, including `summary/tables`.
- `grep -nE '1_raw|2_lncRNA|1_lncRNA|PROJECTION_TOKEN|split_sample|\.from_' bin/run_summary_tables.py`
  returns nothing.
- Results record the runtime of RUN_SUMMARY_TABLES on V2, before and after.

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|
| D27h-1 | A2 | the record index reads the records' values | gffstats values go through the same row builders as the MultiQC tables (`fomo_stats.gffstats_transcript_rows` / `gffstats_gene_rows`, moved there from `gff_stats_to_mqc.py`, which imports them) | the old path formatted each value for the table and parsed it back (`fmt` then `num`): a mean of 100.004 became "100.00", then 100.0, then "100". Reading the raw JSON value would print "100.00" and change `summary/tables` bytes | 27i deletes `gff_stats_to_mqc.py`; the builders stay in `fomo_stats.py` for the renderer (updated by D27h-1) |

## Results

- `run_summary_tables.py` reads only `--records` (a directory of fomo-data/1 records) and the
  roles CSV. `RecordIndex` builds each index once (`cached_property`): annotation rows by
  (species, side, kind, track), with the target-side raw stats falling back to the source side of
  the same species; projected models by (target, source, track); TD2 by (stage, species); the
  transcript funnel; alignment; accuracy by (target, source, track); curation by (target, track).
  Each section is `compute_<section>` (pure) plus an emit step; `RunContext` holds roles, enabled
  feature types and switches; `check_funnel_consistency()` holds the three cross-checks; errors
  are `FomoError`; `fmt` is the shared one in `fomo_stats`.
- `grep -nE '1_raw|2_lncRNA|1_lncRNA|PROJECTION_TOKEN|split_sample|\.from_' bin/run_summary_tables.py`
  is empty; no per-file ignore; no function above complexity 10 (was `section_aggregates` 11
  with C901 ignored); 882 → 738 lines.
- Regression: `tests/run_summary_tables/fixtures/{V1,V2}/records` (every record the V1/V2 gate runs
  wrote: 26 and 74 files) reproduce the golden `expected/` tables byte for byte; the old `mqc/`
  inputs are deleted. RUN_SUMMARY gets the union of the subworkflows' `data` emits.
- RUN_SUMMARY_TABLES on V2: 0.41 s, 15.1 MB before; 0.45 s, 13.9 MB after. MULTIQC_RUN_SUMMARY
  renders every section (checked in the V2 report).
- Gates: G1 278 passed; G2 pass; G3 S1–S4 pass; G4 V1–V4 identical, `summary/tables` included.
  Cluster use 0.3 CPU-h.
