# 27i — MultiQC rendered from records; adapters, duplicates and the old TSV path deleted

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decisions: U5 (duplicate
target stats), U7, T2. Review items: P2, P3, D2, D4, N3, and F8/F16/F25 of the module review.

## Status: done (2026-10-15)

## Goal

All MultiQC custom content is produced in **one place**: `bin/render_mqc.py`, which runs inside the
per-target MULTIQC task (T2), from records and a **section registry** that a test keeps
consistent with the MultiQC YAMLs. Today's adapter processes, their scripts, the `*_mqc.tsv`
outputs of the producers and the three GFF-stats modules collapse. The deliverables and the run
summary stay identical.

## Steps

### A. Section registry

- [x] **A1. Write `bin/fomo_sections.py`.** For each per-target section it records:
  - the record `section` + filter (e.g. `side`, `kind`);
  - the columns (key, title, format, hidden);
  - the row-label function, which turns identity into a display label (the ordering ladders now
    exist **only here**, as presentation);
  - the output filename suffix;
  - the MultiQC section id.

  Sections: gffstats_input, gffstats_genes, gffstats_projection, fasta_stats, samtools_align,
  transcript_filter, td2_coding, curation_observed/discarded/novel, and the accuracy scatter.
- [x] **A2. Write `tests/test_multiqc_registry.py`.** Every registry suffix has an `sp` pattern in
  `assets/multiqc/sections.yml`. Every section id has an entry in `report_section_order`
  (`main.yml`). Every `sp` pattern has a registry entry. The same check covers the run-summary
  pair against `run_summary_tables.py`'s output names.

### B. Renderer

- [x] **B1. Write `bin/render_mqc.py --records DIR --target T --out DIR`.**
  - It renders every registry section from that target's records plus the broadcast
    (source-side) records.
  - It also renders the accuracy scatter. That replaces `gffcompare_accuracy_mqc.py` and
    GFFCOMPARE_TO_SCATTER; the 10-colour palette that wraps for large runs becomes a
    hue-cycled or source-ranked palette.
  - Stdlib only. It must run on the MultiQC image's Python 3.14; pytest runs it on 3.11.
- [x] **B2. Tests.** (updated by D27i-3)
  - Per section: the rendered TSV for the V2 fixtures has the same rows (labels) and values as
    today's adapter TSVs (fixtures from 27g E1).
  - The only allowed difference is the duplicate `_target` rows that U5 removes.
- [x] **B3. Rework the MULTIQC module script** to run `render_mqc.py …` and then `multiqc …`.
  - Inputs: records, native tool files (samtools `.stats` for MultiQC's native module, the
    gffcompare `.stats` for the native gffcompare section), the software-versions YAML (27c),
    and the configs. (updated by D27c-5: MULTIQC must not emit to the versions topic, and the
    YAML must not reach any process that does — e.g. RUN_SUMMARY_TABLES.)
  - MULTIQC_RUN_SUMMARY does not render; its inputs are already custom content.

### C. Delete the old path

- [x] **C1. Delete the modules:** GFF_STATS_TO_MQC (both aliases), SAMTOOLS_TO_MQC, SEQKIT_TO_MQC,
  GFFCOMPARE_TO_SCATTER.
- [x] **C2. Delete the scripts:** `gff_stats_to_mqc.py` (its row builders live in `fomo_stats.py`; updated by D27h-1), `samtools_to_mqc.py`, `seqkit` adapter
  awk, `curation_to_mqc.py`, `gffcompare_accuracy_mqc.py`. Delete `fomo_stats.FNAME_RE` and
  `parse_identity` if nothing else uses them; the accuracy identity now comes from records.
- [x] **C3. Producers stop writing `*_mqc.tsv`:** FILTER_TRANSCRIPT (`--counts`), TD2_CODING_FILTER
  (`--mqc`), CURATE_MODELS (the `curation_to_mqc` call). `*.curation.{tsv,json}` are
  deliverables and stay.
- [x] **C4. Merge the three GFF-stats modules into one, `GFF_STATS`.** (updated by D27i-1)
  - Input `tuple val(meta), path(gffs), path(manifest)`. It loops sequentially, under strict
    mode, over the manifest rows and runs gff-feature-stats plus the printf record wrapper.
    Image `img_ubuntu`, label `amd64_only`.
  - Calls:
    - source side: one task per species (raw + filtered per feature type + decoys per feature
      type);
    - target side: **only pure targets that carry a gff3 (U5)**;
    - projected: one task per (target, track).
  - The manifests are built in the subworkflows from metas: one identity row per file.
- [x] **C5. Remove the config blocks that only existed for the deleted modules:** ext.sample_name
  ladders, role/section knobs, the `.target` segment.
- [x] **C6. Rewrite the MultiQC YAMLs.**
  - `sections.yml`: `custom_data` blocks generated from, or checked against, the registry.
  - `main.yml`: a documented order numbering (hundreds per pipeline stage), entries for the
    accuracy scatter and the native gffcompare section, and `extra_fn_clean_exts` updated.
  - Text in present tense.
- [x] **C7. Simplify REPORTING.**
  - It routes records by `meta.target_id`, broadcasting the source-side ones, and passes
    records, native files and versions to MULTIQC.
  - The "flatten before combine" handling disappears where every item is one file; keep `asList`
    only where a list-valued emit remains. (updated by D27e-1: there is no `asList`; the one
    flattening step is `[p].flatten()` in `routeByTarget` and RUN_SUMMARY.)

### D. Tests and counts

- [x] **D1. Update the stub matrix S1–S3 task counts** (fewer tasks: no adapters, one GFF_STATS
  per species, no target duplicates). Record the old and new counts per process in Results.
- [x] **D2. Run the real V2 run and compare its trace task counts with 27a-A9.** Record them in
  Results; this is the "fewer tasks" goal.

## Acceptance

- G1, G2 (no per-file ignores left anywhere), G3 (with the updated counts), and G4 identical
  (deliverables + `summary/tables`).
- **MultiQC check** (informational under U1, but required here): for V1–V3, every section id in
  the baseline `snapshot.info.tsv` (27a) is still present with the same row count. The only
  exceptions are the U5 removal of duplicate `_target` rows (fewer rows, listed) and the new
  software-versions section.
- No module named `*_TO_MQC` or `GFFCOMPARE_TO_SCATTER` remains. One GFF_STATS module remains.
- `grep -rn '_mqc.tsv' modules` lists only RUN_SUMMARY_TABLES and the renderer outputs inside
  MULTIQC.

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|
| D27i-1 | C4 | source side: one GFF_STATS task per species; manifests of identity rows built in the subworkflows | one module whose loop runs over its rows: a manifest (`file, source`, the projected batches, from 27g) or the task's one GFF3 (source side and pure targets, one task per file, the task's identity) | grouping a species' raw, filtered and decoy GFF3s needs a per-species size key and holds each species' stats until its slowest decoy chain, for a stage that grows only with S (at most S·(1+2F) tasks); per-row identity manifests would need a second definition of the identity next to `ext.identity` | none |
| D27i-2 | B3 | MULTIQC inputs: records, native files, versions, configs | a `records` input slot (staged in `records/`); MULTIQC renders only when it has records, so MULTIQC_RUN_SUMMARY passes `[]` | one module serves both reports | none |
| D27i-3 | B2 | rendered TSVs compared with today's adapter TSVs (fixtures from 27g E1) | `tests/records/test_render_mqc.py` checks all 121 captured table/record pairs (labels in order, values equal); the 27g bridge test and its capture script are deleted with the adapters, the fixtures stay as the oracle | the adapters no longer exist to produce tables, so the captured pairs are the only reference | none |

## Results

> **Correction (found in 27k D4):** C5 left the CURATE_MODELS `saveAs` filter on `_mqc.tsv`, a file
> it no longer writes, and a `fomo_stats.py` comment still named `curation_to_mqc.py`. 27k removed
> both; the published paths do not change (G4).

- **One place.** `bin/fomo_sections.py` (registry: 10 table sections + the accuracy scatter, the
  row labels and their order tokens) and `bin/render_mqc.py` (runs in MULTIQC on Python 3.14,
  stdlib only). The scatter colours the allModels aggregate dark grey and gives every species its
  own evenly spaced hue (distinct for 40 sources, tested). `tests/test_multiqc_registry.py` ties
  the registry to `sections.yml`/`main.yml` and `run_summary_tables.py`'s outputs to the
  run-summary pair.
- **Deleted.** Modules GFF_STATS_TO_MQC (both aliases), SAMTOOLS_TO_MQC, SEQKIT_TO_MQC,
  GFFCOMPARE_TO_SCATTER, GFF_STATS_PROJECTED_BATCH; scripts `gff_stats_to_mqc.py` (its row
  builders live in `fomo_stats.py`), `samtools_to_mqc.py`, `curation_to_mqc.py`,
  `gffcompare_accuracy_mqc.py`; the producers' `*_mqc.tsv` (`--counts/--sample` of
  filter_transcript.py, the TD2 table and `--sample/--feature`, the `curation_to_mqc` call);
  `fomo_stats.FNAME_RE`, `parse_identity`, `Identity`, `feature_class`, `transcript_sn_pr`; every
  per-file lint ignore; the config blocks of the deleted modules, the sample-name ladders and the
  `.target` segment. Target-side stats run only for pure targets with a gff3 (U5).
- **Acceptance greps.** No `*_TO_MQC` or `GFFCOMPARE_TO_SCATTER` module; one GFF_STATS module (three
  aliases); `_mqc.tsv` in `modules/` only in RUN_SUMMARY_TABLES and MULTIQC's renderer call.
- **MultiQC check** (V1–V3 against the 27a `snapshot.info.tsv`): every table keeps its row count
  except (a) the new Software Versions table (V1 23–24 rows per report, 22 in the run summary;
  V2 28–29 / 27) and (b) in alciphron's report (a `both` species) `gffstats_input` 17 → 11 (V1/V3)
  and 24 → 17 (V2), `gffstats_genes` 9 → 6 and 16 → 12: exactly its `_target` rows
  (`1_raw_target.*`, `2_lncRNA_target.*`, `4_mRNA_target.*`), nothing added. Hippothoe (a pure
  target) keeps its target rows.
- **Task counts** (V2 trace, 27a → 27i): 167 → 133. Gone: GFF_STATS_TO_MQC 10, SEQKIT_TO_MQC 8,
  SAMTOOLS_TO_MQC 8, GFF_STATS_TARGET_TO_MQC 4, GFFCOMPARE_TO_SCATTER 1; GFF_STATS_TARGET 4 → 1;
  GFF_STATS_PROJECTED_BATCH 8 → GFF_STATS_PROJECTED 8. Stub matrix: S1/S2 GFF_STATS_TARGET 3/4 → 1,
  S3 2 → 0, the five adapter processes removed.
- **Gates:** G1 266 passed; G2 pass (no per-file ignores anywhere); G3 S1–S4 pass with the new
  counts; G4 V1–V4 identical (deliverables and `summary/tables`). Cluster use 0.3 CPU-h.
