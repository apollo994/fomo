# 15b — Phase B: raw and collapsed aggregates

Sub-plan of [`15_single_mapping.md`](15_single_mapping.md). Requires Phase A
([`15a`](15a_single_mapping_align.md)) checkpoint to have passed.

## Goal

Stop building the consensus by handing `gffcompare` N per-source files. Build it from the
single `allModels.raw.gff3` instead, add `top3.raw.gff3` as a subset of the same file, and
score all four aggregates. After this phase an annotated target yields **S + 4** accuracy
values: per-source, `allModels.raw`, `allModels.collapsed`, `top3.raw`, `top3.collapsed`.

Out of scope: publish layout, MultiQC closures, resources, docs — all Phase C.

## Changes

### `subworkflows/local/projection.nf` — collapse from `allModels`

Replace the `ch_projected_persource.filter { meta.id != meta.target_id }.map{…}.groupTuple()`
block feeding `GFFCOMPARE_COMBINE` with the single filtered GFF:

```groovy
GFFCOMPARE_COMBINE(
    ch_allmodels_raw.map { meta, gff -> tuple(meta + [gtype: meta.feature_type + (meta.decoy ? '_decoy' : '')], gff) },
    [[:], [], []],   // no reference sequence  (-s)
    [[:], []]        // no reference annotation (-r) → pure combine mode
)
COMBINED_GTF_TO_GFF(
    GFFCOMPARE_COMBINE.out.combined_gtf.map { meta, gtf -> tuple(meta + [id: 'allModels_collapsed'], gtf) },
    []
)
```

**Self is now included** (plan 15 decision) — the `meta.id != meta.target_id` filter goes away
here. The self models are in `allModels.raw`, so they are in the collapse. Delete the
"self-pairs are EXCLUDED here" comment block and replace it with a note saying self is kept in
`allModels` but still excluded from the top-N ranking pool, so the two rules are not confused
later.

Emit shape becomes:

```groovy
emit:
bam                 = MINIMAP2_ALIGN.out.bam
index               = MINIMAP2_ALIGN.out.index
allmodels_raw       = ch_allmodels_raw                       // id = 'allModels_raw'
allmodels_collapsed = COMBINED_GTF_TO_GFF.out.gffread_gff    // id = 'allModels_collapsed'
gff3                = ch_projected_persource
                        .mix(ch_allmodels_raw.map { m, g -> tuple(m + [id: 'allModels_raw'], g) })
                        .mix(COMBINED_GTF_TO_GFF.out.gffread_gff)
mqc_files           = ch_mqc_files
```

`gff3` is what `BENCHMARKING` consumes, so `allModels.raw` and `allModels.collapsed` get
scored with **no change to `benchmarking.nf`** — its `[target_id, feature_type]` inner join
already handles them, and drops them for targets with no `gff3`.

`GFF_STATS_PROJECTED` should cover the aggregates too — run it on the union, not just the
per-source channel.

### `subworkflows/local/consensus_top.nf` — subset, then collapse

`SELECT_TOP_SOURCES` and its per-target grouping stay **exactly** as they are, including the
`meta.id != meta.target_id` self-exclusion and the `!(meta.id in ['combined','top3'])` guard
(extend that list to the new aggregate ids). Rationale in plan 14 still holds: self always
wins on F1, and `bin/select_top_sources.py` keys its score map by source alone, so the task
must stay per-target.

Replace the `ch_projected_real ⋈ ch_selected → groupTuple → GFFCOMPARE_COMBINE_TOP` block:

```groovy
ch_subset_in = ch_allmodels_raw                       // new take: input, real tracks only
    .filter { meta, _gff -> !meta.decoy }
    .map { meta, gff -> tuple(meta.target_id, meta, gff) }
    .combine(SELECT_TOP_SOURCES.out.csv.map { meta, csv -> tuple(meta.target_id, csv) }, by: 0)
    .map { _tid, meta, gff, csv -> tuple(meta + [id: 'top3_raw'], gff, csv) }

SUBSET_GFF_BY_SOURCE(ch_subset_in)                    // alias of GFF_BY_SOURCE, --keep mode
ch_top_raw = SUBSET_GFF_BY_SOURCE.out.gff3.transpose()

GFFCOMPARE_COMBINE_TOP(ch_top_raw.map { m, g -> tuple(m + [gtype: m.feature_type], g) }, [[:], [], []], [[:], []])
COMBINED_TOP_GTF_TO_GFF(GFFCOMPARE_COMBINE_TOP.out.combined_gtf.map { m, g -> tuple(m + [id: 'top3_collapsed'], g) }, [])
```

`GFF_BY_SOURCE` needs a third input slot for the keep-list, or the module takes
`tuple val(meta), path(gff3), path(keep)` with `keep` optional (`[]` in split mode). Pick one
shape and use it in both aliases.

Score **both** aggregates — today only the collapsed model is benchmarked:

```groovy
ch_top_scored = ch_top_raw.mix(COMBINED_TOP_GTF_TO_GFF.out.gffread_gff)
GFF_STATS_TOP(ch_top_scored.map { m, g -> tuple(m + [kind: 'projected'], g) })
// ⋈ target reference on [target_id, feature_type] — NEVER feature_type alone
GFFCOMPARE_TOP(...)
```

The projection ⋈ reference join must keep the `[target_id, feature_type]` key, for the reason
`CLAUDE.md` documents: on `feature_type` alone every projection is also compared against every
*other* target's annotation, and since `ext.prefix` is built from the query meta all T tasks
emit the same filename into the same directory.

### `workflows/fomo.nf`

`CONSENSUS_TOP` now takes `PROJECTION.out.allmodels_raw` instead of `ch_projected_real`.
Update the `filter` on lines 81–82 to exclude the new aggregate ids
(`allModels_raw`, `allModels_collapsed`) where `'combined'` was excluded.

### `conf/modules.config`

| Process | `ext.prefix` |
|---|---|
| `FOMO:PROJECTION:GFFCOMPARE_COMBINE` | `${meta.target_id}.allModels.${meta.gtype}` |
| `FOMO:PROJECTION:COMBINED_GTF_TO_GFF` | `${meta.target_id}.allModels.${meta.feature_type}${decoy}.collapsed` |
| `FOMO:CONSENSUS_TOP:SUBSET_GFF_BY_SOURCE` | `ext.args` → `--keep <csv> --name ${meta.target_id}.top3.${meta.feature_type}.raw.gff3` |
| `FOMO:CONSENSUS_TOP:GFFCOMPARE_COMBINE_TOP` | `${meta.target_id}.top3.${meta.gtype}` |
| `FOMO:CONSENSUS_TOP:COMBINED_TOP_GTF_TO_GFF` | `${meta.target_id}.top3.${meta.feature_type}.collapsed` |
| `FOMO:CONSENSUS_TOP:GFF_STATS_TOP$` | `${meta.target_id}.from_${meta.id}.${meta.feature_type}.projected` |

Benchmark `GFFCOMPARE` / `GFFCOMPARE_TOP` keep the
`${target}.from_${meta.id}.${ft}[.decoy].gffcompare` shape. That filename is the contract
parsed by `bin/select_top_sources.py` (`FNAME_RE`) and `bin/gffcompare_accuracy_mqc.py`
(same regex), so with `meta.id ∈ {allModels_raw, allModels_collapsed, top3_raw,
top3_collapsed}` both scripts keep working and the aggregates get their own colours on the
accuracy scatter. Do **not** invent a different naming shape for them.

### `bin/select_top_sources.py`

Extend the guard on line 81 from `source in ("combined", "top3")` to also skip
`allModels_raw`, `allModels_collapsed`, `top3_raw`, `top3_collapsed`. Second layer only — the
channel filter in `consensus_top.nf` is authoritative.

## Checkpoint B

Run `nextflow run main.nf -profile test,crg --outdir results_15b`.

1. **Four aggregate stats files** exist for `Lycaena_alciphron_282377`:
   `*.from_allModels_raw.lnc_RNA.gffcompare.stats`, `*.from_allModels_collapsed.*`,
   `*.from_top3_raw.*`, `*.from_top3_collapsed.*` — plus the 2 per-source ones. Six total.
2. **Four aggregate GFF3s** exist: `allModels.lnc_RNA.raw.gff3`,
   `allModels.lnc_RNA.collapsed.gff3`, `top3.lnc_RNA.raw.gff3`, `top3.lnc_RNA.collapsed.gff3`.
3. **Monotonicity sanity**: transcript count `allModels.collapsed ≤ allModels.raw`, and
   `top3.raw ⊆ allModels.raw` (every `top3.raw` transcript id is present in `allModels.raw`).
4. **Self inclusion is visible**: `allModels.raw.gff3` contains
   `source=Lycaena_alciphron_282377` rows for that target, and `top3.raw.gff3` does **not**
   (self excluded from the ranking pool). `top_sources.csv` must not list the target itself.
5. **Per-source outputs unchanged from Phase A** — `diff -r` the per-source GFF3s and their
   `.stats` between `results_15a` and `results_15b`; nothing should have moved. Phase B only
   adds aggregates.
6. **Un-annotated target** `Lycaena_hippothoe_580924`: `allModels.raw` and
   `allModels.collapsed` present, no `gffcompare/`, no `select_top_sources/`, no `top3.*`.
7. **`select_top_sources.py` did not rank an aggregate** — read
   `select_top_sources/*.top_sources.csv` and confirm only real species names appear.
8. **Re-read against the plans**: confirm the `[target_id, feature_type]` join key survived in
   both `benchmarking.nf` and `consensus_top.nf`; that no subworkflow got nested; that
   `SELECT_TOP_SOURCES` is still one task per target; and that the decision table in plan 15
   (self in `allModels`, self out of ranking) matches what the code does.

## Checkpoint B — RESULT (run 2026-08-12, `-profile test,crg`, `results_15b`)

**PASSED.** 92 tasks (88 cached on a `-resume`; see the note on the OOM below).

| Check | Result |
|---|---|
| 1. Six `.stats` for the annotated target | `from_Lycaena_alciphron_282377`, `from_Lycaena_thersamon_265360`, `from_allModels_raw`, `from_allModels_collapsed`, `from_top3_raw`, `from_top3_collapsed` ✓ |
| 2. Four aggregate GFF3s | `allModels.lnc_RNA.raw`, `allModels.lnc_RNA.collapsed`, `top3.lnc_RNA.raw`, `top3.lnc_RNA.collapsed` ✓ |
| 3. Monotonicity | allModels 124 → collapsed **105**; top3.raw 64 → collapsed **63**; `top3.raw ⊆ allModels.raw` (0 ids outside) ✓ |
| 4. Self | `allModels.raw` sources = {alciphron, thersamon}; `top3.raw` sources = {thersamon}; `top_sources.csv` lists only thersamon ✓ |
| 5. Per-source unchanged vs `results_15a` | 4/4 GFF3s byte-identical, 2/2 per-source `.stats` content-identical ✓ |
| 6. Un-annotated target | `filter_allmodels`, `combined_gtf_to_gff`, `split_gff_by_source`, `minimap2_align`, `multiqc` only — no `gffcompare/`, no `select_top_sources/`, no `top3.*`; no `targets/null` ✓ |
| 7. No aggregate ranked | `top_sources.csv` = one real species ✓ |
| 8. Plan/contract re-read | `by: [0, 1]` intact in `benchmarking.nf:52` and `consensus_top.nf:105`; selectors still two levels deep; `SELECT_TOP_SOURCES` = 1 task; `meta.id != meta.target_id` now appears in exactly one place (the ranking filter) ✓ |

Task counts: `GFFCOMPARE_COMBINE` 2, `SUBSET_GFF_BY_SOURCE` 1, `GFFCOMPARE_COMBINE_TOP` 1,
`GFFCOMPARE_TOP` **2** (was 1), `BENCHMARKING:GFFCOMPARE` **4** (was 3),
`GFF_STATS_PROJECTED` **8** (was 6).

### The four aggregates, transcript-level Sn / Pr

| Model | Sn | Pr |
|---|---|---|
| `Lycaena_alciphron_282377` (self, ceiling control) | 88.1 | 98.3 |
| `Lycaena_thersamon_265360` | 26.9 | 28.1 |
| `allModels_raw` | 88.1 | 47.6 |
| `allModels_collapsed` | 88.1 | **56.2** |
| `top3_raw` | 26.9 | 28.1 |
| `top3_collapsed` | 26.9 | **28.6** |

Collapsing costs no sensitivity and buys precision (47.6 → 56.2 on allModels) — exactly the
raw-vs-collapsed comparison this phase exists to produce. All six land on the accuracy scatter
with distinct colours, and all six appear in the projection GFF-stats table.

Two artefacts of the tiny test set, not of the code: `top3_raw` is identical to
`Lycaena_thersamon_265360` because with 2 sources and self excluded the top-N *is* thersamon
alone; and `allModels`'s Sn jumps to the self ceiling now that self is included (plan 15
decision), which is why its numbers are not comparable to the old `from_combined` (26.9/28.6).

### Notes

- **The run was OOM-killed at ~98% and resumed.** Not a pipeline fault: the interactive session
  job is capped at 2 GB and the Nextflow JVM exceeded the cgroup (exit 137) after submitting the
  final reporting tasks. Re-running with `NXF_OPTS='-Xms256m -Xmx1200m' -resume` finished in
  1m37s (4 new tasks, 88 cached). Use that `NXF_OPTS` for any run driven from a 2 GB session.
- **`Skipping output binding … optional files are missing` is pre-existing and benign.** It fires
  37 times across `EXTRACT_SEQUENCES`, `MINIMAP2_ALIGN`, `MULTIQC` and every `GFFCOMPARE*` alias,
  because those modules declare several mutually exclusive optional outputs (a `-r` comparison
  never writes `.combined.gtf`; a combine run never writes `.annotated.gtf`). All 6 expected
  files exist for both `top3_raw` and `top3_collapsed`.
- **For Phase C:** `bin/gffcompare_accuracy_mqc.py` still special-cases the literal source
  `"combined"` for stable colour ordering (lines 122–125). That id no longer exists, so the
  branch is dead and colours are now purely alphabetical. Point it at the aggregate ids instead.
