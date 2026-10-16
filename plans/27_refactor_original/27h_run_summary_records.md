# 27h — The run summary reads data records

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decision: U7. Review items:
P2, P4, D2, R5, and F5/F6/F8–F13/F16 of the Python review.

## Status: not started

## Goal

`RUN_SUMMARY_TABLES` builds every table and `run_summary.json` from **records only** (27g) plus the
roles CSV. It no longer parses Sample labels (`1_raw`, `2_lncRNA`, `PROJECTION_TOKEN`),
filename suffixes or `.from_` grammars. The code is split into cached record indexes, pure
`compute_<section>` functions and `emit_<section>` writers. **Outputs stay byte-identical (U1).**

## Steps

- [ ] **A1. Extend the regression fixture from 27a** (`tests/run_summary_tables/fixtures/{V1,V2}/`)
  with the records of the same runs (captured in 27g E1). The existing test (old inputs →
  golden outputs) keeps passing until A5.
- [ ] **A2. Write the new input layer, `RecordIndex`.**
  - It loads every `*.data.json` once and builds the indexes: annotation stats by
    (species, side, kind, track); transcript funnel; TD2; alignment; accuracy by
    (target, source, track); curation by (target, track).
  - Use `functools.cached_property`, so no index is rebuilt; today `accuracy()` is rebuilt 3×
    and `annotation_index()` 3×.
  - **Target-side raw stats:** use the `side=target` record if present, otherwise the
    `side=source` record of the same species (identical file). This prepares 27i's removal of
    duplicate target stats (U5) without changing values.
- [ ] **A3. Split every section into `compute_<section>(ctx) -> records` (pure) and
  `emit_<section>(records)`.**
  - Sections: overview, species, source_funnel, target_types, heatmap, aggregates, donors.
  - The funnel consistency checks become `check_funnel_consistency()`.
  - A `RunContext` holds sources, targets and enabled feature types, derived once.
  - Remove the 4 unused `args` parameters and the out-parameter `summary` dict mutation; `main`
    assembles the JSON.
- [ ] **A4. Errors and naming.**
  - Raise `FomoError` (from `fomo_io`), not `sys.exit` inside the class; one policy for missing
    inputs.
  - Use the shared `fmt`/`pct` helpers and the shared NA policy.
  - Delete `PROJECTION_TOKEN`, `split_sample`, `_curation_key` and every hard-coded `"1_raw"`,
    `"2_lncRNA"` or `"1_lncRNA"`.
- [ ] **A5. Switch `RUN_SUMMARY_TABLES`'s input.**
  - `run_summary.nf` passes the records (the union of the `data` emits) instead of `ch_all_mqc`.
  - Delete the old-input regression test, and rename the record-based one to be the regression
    test. Its expected outputs are the same golden `summary/tables`.
- [ ] **A6. Remove the per-file ignores of `run_summary_tables.py`.** C901 ≤ 15 holds.
- [ ] **A7. MultiQC_RUN_SUMMARY still renders.** Its inputs (`run_*_mqc.*`) are unchanged in name
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

## Results

*(fill in)*
