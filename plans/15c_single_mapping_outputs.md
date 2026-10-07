# 15c — Phase C: output layout, reporting, resources, docs

Sub-plan of [`15_single_mapping.md`](15_single_mapping.md). Requires Phase B
([`15b`](15b_single_mapping_aggregates.md)) checkpoint to have passed.

## Goal

The pipeline is functionally correct after Phase B but its output tree is still named after
processes that no longer mean anything (`filter_allmodels/`, `combined_gtf_to_gff/`), MultiQC
row labels still say `from_allModels` where the merged BAM is concerned, and `MINIMAP2_ALIGN`
is sized for one source rather than all of them. Phase C makes the deliverable presentable and
the docs true.

## Changes

### Output layout — `nextflow.config` + `conf/modules.config`

The publish block in `nextflow.config` (lines 59–67) derives the directory from the process
name, which was fine when `BAM_TO_GFF` / `COMBINED_GTF_TO_GFF` *were* the deliverables. Replace
it with explicit per-process `publishDir` entries in `conf/modules.config`, keeping the
`targets/<target>/` root and the `saveAs` that drops `versions.yml`:

```
targets/<target>/alignment/    <target>.allModels.<ft>[.decoy].raw.bam
                               <target>.allModels.<ft>[.decoy].raw.bam.bai
targets/<target>/annotation/   <target>.allModels.<ft>[.decoy].raw.gff3
                               <target>.allModels.<ft>[.decoy].collapsed.gff3
                               <target>.top3.<ft>.raw.gff3
                               <target>.top3.<ft>.collapsed.gff3
                               <target>.from_<source>.<ft>[.decoy].projected.gff3
targets/<target>/gffcompare/   *.stats *.tracking *.loci *.tmap *.refmap
targets/<target>/select_top_sources/   <target>.top_sources.csv
targets/<target>/multiqc/      multiqc_report.html + multiqc_report_data/
```

This is exactly what the prompt asks for — "all gff annotations and the allModels.raw.bam".
Keep the `meta.target_id` closure: it is still the per-target segment, and a process that
loses the key publishes to `targets/null/`.

Retain the sizing note from the old block, updated: at 20 targets × F × D this now publishes
~40 BAMs instead of ~1600, so `MINIMAP2_ALIGN` no longer needs a disk caveat.

### MultiQC — `conf/modules.config` closures

`SAMTOOLS_STATS` / `SAMTOOLS_TO_MQC` now describe the merged alignment, one row per
(target, class). Update `ext.sample_name` to `${meta.target_id}.allModels.<order>` using the
same sparse-safe `1_lncRNA / 2_decoy_lncRNA / 3_mRNA / 4_decoy_mRNA` ladder — a disabled class
simply produces no row, and the `'unknown'` fallback stays unreachable.

Add the new filename suffixes to `extra_fn_clean_exts` in `assets/multiqc/main.yml`
(`.allModels.lnc_RNA.raw`, `.collapsed`, `.top3.lnc_RNA.raw`, …) so sample names render
readably. No new `custom_data` block or `sp` pattern is needed — the aggregates reuse the
existing `gffstats_projection` and gffcompare sections. Do **not** add an `# id:` header to
any adapter output, and do not combine `custom_content.order` with `report_section_order`;
either mistake double-renders a section.

Document the accepted regression in the report comment: alignment metrics are per-target-class,
not per-source; per-source model counts remain in the projection GFF-stats table.

### Resources — `conf/base.config` / `conf/crg.config`

`MINIMAP2_ALIGN` now aligns the union of every source in one task. Memory is dominated by the
target index rather than the read set, so `process_medium`'s 24 GB should still hold, but wall
time scales with total read count. Give `MINIMAP2_ALIGN` its own `withName` block with a
raised `time` rather than moving the whole `process_medium` label — other processes share it.

Check the `-profile test` `resourceLimits` (`[cpus: 2, memory: '4.GB', time: '30.min']`) still
admits the merged alignment on the subsampled test data.

### Documentation

- **`CLAUDE.md`** — update the task-count formulas in *Feature tracks* (`F·D·T` alignments and
  BAM→GFF, not `F·S·D·T`; add `MERGE_SOURCE_FASTA` at `F·D`; benchmark count becomes
  `(F·S·D + 2·F·D)·T_g`); update the *Subworkflow map* table (PROJECTION now merges then
  splits; CONSENSUS_TOP subsets rather than group-combines); rewrite the self-pair paragraph —
  self is **kept** in `allModels` and still **excluded** from the top-N ranking pool; note the
  new `source=`-attribute split contract and the duplicate-ID guard; update the results layout
  paragraph to `alignment/` + `annotation/`.
- **`plans/15_single_mapping.md`** — mark the phases done.
- **`README.md`** if it documents the output tree.

## Checkpoint C

Run `nextflow run main.nf -profile test,crg --outdir results_15c`, then with
`--include_mrna --include_decoy` into `results_15c_full`.

1. **Layout matches the table above** exactly, for both targets. No `targets/null/`. No
   leftover `bam_to_gff/`, `combined_gtf_to_gff/`, `combined_top_gtf_to_gff/`,
   `filter_allmodels/` directories.
2. **Content unchanged from Phase B** — the files moved, they did not change:
   `diff` each `results_15b` GFF3/stats against its new location; all must be identical.
3. **One MultiQC report per target**, both render. The accuracy scatter shows the per-source
   dots **plus** four aggregate dots with distinct colours and readable tooltips. No section
   appears twice. The samtools table has one row per (target, class), correctly labelled.
4. **`--include_mrna --include_decoy` run**: 4 merged FASTAs, `4 × T` alignments, an
   `allModels` + `top3` set per real class, decoy classes get `allModels` but no `top3`
   (decoys are filtered out of the top-N path). `mRNA` and decoy per-source GFF3s still
   byte-identical to the Phase A baseline.
5. **Resources**: no `MINIMAP2_ALIGN` retry in `execution_trace.txt` for either run.
6. **Params/schema in sync**: `nextflow run . --help` renders, no `failUnrecognisedParams`
   abort. Any param added in Phases A–C is in both `nextflow.config` and
   `nextflow_schema.json`.
7. **Re-read against the plans**: walk `CLAUDE.md` against the code — every formula, every
   meta key in the contract (`feature_type`, `decoy`, `target_id`, `kind`), the "never nest
   the subworkflows" rule, and the `[target_id, feature_type]` join rule. Then walk plan 15's
   decision table and cardinality table against `execution_trace.txt`. Anything that no longer
   matches is either a doc bug to fix or a code bug to fix — do not leave it as a discrepancy.
8. **Accept and record the two known deltas** from plan 15: per-source lncRNA GFF3s differ
   from `results_baseline/` by TD2-dropped transcripts only, and per-source alignment rows are
   gone from MultiQC. Confirm both are still exactly those two, and nothing more.

## Checkpoint C — RESULT (run 2026-08-12)

**PASSED, all 8 items.** `results_15c` (92 tasks, 13m07s) and `results_15c_full`
(`--include_mrna --include_decoy`, 263 tasks, 15m37s).

| Check | Result |
|---|---|
| 1. Layout | Exact. `alignment/` 2, `annotation/` 6 (annotated) / 4 (un-annotated), `gffcompare/` 36, `select_top_sources/` 1, `multiqc/`. No `bam_to_gff/`, `combined_gtf_to_gff/`, `combined_top_gtf_to_gff/`, `filter_allmodels/`, `split_gff_by_source/`, `subset_gff_by_source/`, `minimap2_align/`. No `targets/null/` |
| 2. Content unchanged from Phase B | **51/51 files byte-identical** — moved, not changed |
| 3. MultiQC | 2 reports; `.raw` stripped (`<target>.allModels.lnc_RNA`); 6 scatter points, aggregates on fixed palette slots 1–4; every section anchor appears exactly **once** |
| 4. Both flags | `MERGE_SOURCE_FASTA` **4** (F·D); `MINIMAP2_ALIGN`/`BAM_TO_GFF`/`FILTER_ALLMODELS`/`SPLIT_GFF_BY_SOURCE`/`GFFCOMPARE_COMBINE` **8** each (F·D·T); `TD2_PREDICT` **2** (T); bench `GFFCOMPARE` **16** ((F·S·D+2·F·D)·T_g); `SUBSET_GFF_BY_SOURCE` **2** (F·T_g); `GFFCOMPARE_TOP` **4** (2·F·T_g). 4 top3 files, **0** decoy top3. All **16** per-source GFF3s byte-identical to the old-code flags baseline |
| 5. Resources | Zero retries in either run |
| 6. Params/schema | `--help` renders, no errors; no params added in A–C so the two files stay in sync |
| 7. Docs vs code | All 13 cardinality formulas match `execution_trace.txt`; all five meta keys (`feature_type` 84, `decoy` 49, `target_id` 51, `kind` 17, `gtype` 3 refs) in active use; `by: [0, 1]` intact in `benchmarking.nf:52` and `consensus_top.nf:105`; every selector still two levels deep |
| 8. The two deltas | Confirmed, and nothing more: per-source lncRNA GFF3s are byte-identical to `results_baseline` (TD2 dropped nothing on this data — see 15a); samtools rows 2 → **1**, while projection GFF-stats rows went 4 → **6** |

### Changes made

- **`nextflow.config`** — the process-name-derived publishDir alternation is gone, replaced by
  a pointer comment plus the two traps worth knowing (modules.config wins by include order;
  `enabled` only honours a static boolean).
- **`conf/modules.config`** — new "Result publishing" section at the bottom, grouped by
  destination. Five processes feed `annotation/`; their names are disjoint by construction.
- **`conf/base.config`** — `MINIMAP2_ALIGN` gets `time = { 1.h * task.attempt }`, sized from
  measurement: longest per-source alignment 10.5 s vs longest merged 21.2 s at S = 2, so
  roughly linear in S. Memory deliberately **not** raised — peak RSS moved only 2.5 GB → 2.6 GB
  (both on the mRNA track), confirming the target index dominates, not the read set. `time`
  only, because minimap2 interpolates `task.cpus` into its command line.
- **`assets/multiqc/main.yml`** — `.raw` added to `extra_fn_clean_exts`; `report_comment` now
  states that Samtools describes one merged alignment and explains the raw/collapsed
  aggregates.
- **`bin/gffcompare_accuracy_mqc.py`** — the dead `"combined"` colour special-case is replaced
  by the four aggregate ids, so the headline models keep stable colours across runs.
- **`CLAUDE.md`** — cardinality formulas, subworkflow map, self-pair rule, meta contract
  (`gtype`, the pseudo-source `id` values, the three places the aggregate-id list must agree),
  output tree, the TD2-pooling note, a new **merge/split contract** section, and the
  per-target-class alignment-metrics note.

### Deviations from this plan

- **`extra_fn_clean_exts` gained only `.raw`, not the full list of new suffixes.** Every other
  candidate (`.collapsed`, `.projected`, `.top3.lnc_RNA.raw`) would have been dead config: those
  files reach MultiQC through `*_TO_MQC` adapters that set the sample name explicitly, so
  filename cleaning never touches them. `.raw` is the only one that changes a rendered label
  (the native Samtools sections). Verified by reading every label in
  `multiqc_report_data/` before and after.
- **The `MINIMAP2_ALIGN` time raise is precautionary, not demonstrated.** Nothing retried in
  any run, and `-profile test` clamps `time` to 30 min anyway, so the 1 h only takes effect on
  a real clade.
