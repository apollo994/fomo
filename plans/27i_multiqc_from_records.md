# 27i — MultiQC rendered from records; adapters, duplicates and the old TSV path deleted

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decisions: U5 (duplicate
target stats), U7, T2. Review items: P2, P3, D2, D4, N3, and F8/F16/F25 of the module review.

## Status: not started

## Goal

All MultiQC custom content is produced in **one place**: `bin/render_mqc.py`, which runs inside the
per-target MULTIQC task (T2), from records and a **section registry** that a test keeps
consistent with the MultiQC YAMLs. Today's adapter processes, their scripts, the `*_mqc.tsv`
outputs of the producers and the three GFF-stats modules collapse. The deliverables and the run
summary stay identical.

## Steps

### A. Section registry

- [ ] **A1. Write `bin/fomo_sections.py`.** For each per-target section it records:
  - the record `section` + filter (e.g. `side`, `kind`);
  - the columns (key, title, format, hidden);
  - the row-label function, which turns identity into a display label (the ordering ladders now
    exist **only here**, as presentation);
  - the output filename suffix;
  - the MultiQC section id.

  Sections: gffstats_input, gffstats_genes, gffstats_projection, fasta_stats, samtools_align,
  transcript_filter, td2_coding, curation_observed/discarded/novel, and the accuracy scatter.
- [ ] **A2. Write `tests/test_multiqc_registry.py`.** Every registry suffix has an `sp` pattern in
  `assets/multiqc/sections.yml`. Every section id has an entry in `report_section_order`
  (`main.yml`). Every `sp` pattern has a registry entry. The same check covers the run-summary
  pair against `run_summary_tables.py`'s output names.

### B. Renderer

- [ ] **B1. Write `bin/render_mqc.py --records DIR --target T --out DIR`.**
  - It renders every registry section from that target's records plus the broadcast
    (source-side) records.
  - It also renders the accuracy scatter. That replaces `gffcompare_accuracy_mqc.py` and
    GFFCOMPARE_TO_SCATTER; the 10-colour palette that wraps for large runs becomes a
    hue-cycled or source-ranked palette.
  - Stdlib only. It must run on the MultiQC image's Python 3.14; pytest runs it on 3.11.
- [ ] **B2. Tests.**
  - Per section: the rendered TSV for the V2 fixtures has the same rows (labels) and values as
    today's adapter TSVs (fixtures from 27g E1).
  - The only allowed difference is the duplicate `_target` rows that U5 removes.
- [ ] **B3. Rework the MULTIQC module script** to run `render_mqc.py …` and then `multiqc …`.
  - Inputs: records, native tool files (samtools `.stats` for MultiQC's native module, the
    gffcompare `.stats` for the native gffcompare section), the software-versions YAML (27c),
    and the configs. (updated by D27c-5: MULTIQC must not emit to the versions topic, and the
    YAML must not reach any process that does — e.g. RUN_SUMMARY_TABLES.)
  - MULTIQC_RUN_SUMMARY does not render; its inputs are already custom content.

### C. Delete the old path

- [ ] **C1. Delete the modules:** GFF_STATS_TO_MQC (both aliases), SAMTOOLS_TO_MQC, SEQKIT_TO_MQC,
  GFFCOMPARE_TO_SCATTER.
- [ ] **C2. Delete the scripts:** `gff_stats_to_mqc.py`, `samtools_to_mqc.py`, `seqkit` adapter
  awk, `curation_to_mqc.py`, `gffcompare_accuracy_mqc.py`. Delete `fomo_stats.FNAME_RE` and
  `parse_identity` if nothing else uses them; the accuracy identity now comes from records.
- [ ] **C3. Producers stop writing `*_mqc.tsv`:** FILTER_TRANSCRIPT (`--counts`), TD2_CODING_FILTER
  (`--mqc`), CURATE_MODELS (the `curation_to_mqc` call). `*.curation.{tsv,json}` are
  deliverables and stay.
- [ ] **C4. Merge the three GFF-stats modules into one, `GFF_STATS`.**
  - Input `tuple val(meta), path(gffs), path(manifest)`. It loops sequentially, under strict
    mode, over the manifest rows and runs gff-feature-stats plus the printf record wrapper.
    Image `img_ubuntu`, label `amd64_only`.
  - Calls:
    - source side: one task per species (raw + filtered per feature type + decoys per feature
      type);
    - target side: **only pure targets that carry a gff3 (U5)**;
    - projected: one task per (target, track).
  - The manifests are built in the subworkflows from metas: one identity row per file.
- [ ] **C5. Remove the config blocks that only existed for the deleted modules:** ext.sample_name
  ladders, role/section knobs, the `.target` segment.
- [ ] **C6. Rewrite the MultiQC YAMLs.**
  - `sections.yml`: `custom_data` blocks generated from, or checked against, the registry.
  - `main.yml`: a documented order numbering (hundreds per pipeline stage), entries for the
    accuracy scatter and the native gffcompare section, and `extra_fn_clean_exts` updated.
  - Text in present tense.
- [ ] **C7. Simplify REPORTING.**
  - It routes records by `meta.target_id`, broadcasting the source-side ones, and passes
    records, native files and versions to MULTIQC.
  - The "flatten before combine" handling disappears where every item is one file; keep `asList`
    only where a list-valued emit remains.

### D. Tests and counts

- [ ] **D1. Update the stub matrix S1–S3 task counts** (fewer tasks: no adapters, one GFF_STATS
  per species, no target duplicates). Record the old and new counts per process in Results.
- [ ] **D2. Run the real V2 run and compare its trace task counts with 27a-A9.** Record them in
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

## Results

*(fill in)*
