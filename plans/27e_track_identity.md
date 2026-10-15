# 27e — Track registry, meta identity, channel simplification

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decision: T4. Review items:
N1, N2, P7, P10, and F3/F5/F10/F11/F23/F25 of the workflow review.

## Status: done (2026-10-15)

## Goal

A track — feature type plus decoy — is defined **once** and carried in meta. No expression
recomputes it, no task is tagged `allModels`, and channel pairing uses one helper instead of nine
hand-written key/strip patterns.

## Steps

### A. Registry and meta keys

- [x] **A1. Write `subworkflows/local/fomo_utils.nf`** (updated by D27e-1) (it may already exist from 27c-D2) with:
  - `trackOf(String featureType, boolean decoy)` → `[track:, label:, order:, td2:]`:

    | featureType | decoy | track | label | order | td2 |
    |:---|:---|:---|:---|---:|:---|
    | lnc_RNA | false | `lnc_RNA` | `lncRNA` | 1 | true |
    | lnc_RNA | true | `lnc_RNA.decoy` | `decoy_lncRNA` | 2 | false |
    | mRNA | false | `mRNA` | `mRNA` | 3 | false |
    | mRNA | true | `mRNA.decoy` | `decoy_mRNA` | 4 | false |

    An unknown feature type is an error.
  - `withTrack(meta)` → `meta + trackOf(meta.feature_type, meta.decoy)`.
  - `pairById(left, right)`, used in step B.
  - `asList(x)`, which replaces the `instanceof List` idiom.
- [x] **A2. Stamp the track keys once.**
  - PREPROCESSING: where `feature_type`/`decoy` are first set, for the real track
    (`preprocessing.nf` ≈ line 50) and for the decoy track (≈ line 126).
  - PROJECTION's merge: replace `gtype` with the carried keys.
  - Delete `meta.gtype` (underscore form). MERGE_SOURCE_FASTA's prefix becomes
    `allsources.${meta.track}`; that file is intermediate, so U1 is not affected.
- [x] **A3. Replace every `${meta.feature_type}${meta.decoy ? '.decoy' : ''}` with
  `${meta.track}`** (35 sites: `modules.config`, modules, subworkflows). Check with
  `grep -rn "meta.decoy ? '" workflows subworkflows modules conf`, which must return nothing.
- [x] **A4. TD2 branching uses `meta.td2`** (`preprocessing.nf`, `projection.nf`), not
  `!meta.decoy && meta.feature_type == 'lnc_RNA'`.
- [x] **A5. The sample-name ladders in `modules.config` use `meta.order`/`meta.label`** for
  tracks. The kind ladder for raw/filtered stays until 27i deletes these closures. **Sample
  names must stay byte-identical**, because `run_summary_tables.py` parses them until 27h:
  G4 + the run-summary fixture test.

### B. `meta.id` and tags

- [x] **B1. Stop setting `meta.id` to pseudo values (`'allModels'`, `'curated'`).**
  - Target-scoped tasks carry `id = target_id`.
  - Every prefix that relied on the pseudo id is spelled out with `meta.target_id`/`meta.track`
    in `modules.config` (most already are).
  - The pseudo-source `allModels_raw` in the stats grammar is untouched; 27g replaces it.
- [x] **B2. One tag rule** (master §8) in every module. Check: in the V2 trace, no two tasks of
  the same process share a tag. MINIMAP2_ALIGN shows `<target>.<track>`. (updated by D27c-7: 27c
  left every module on its existing tag form, so this step covers all 27 modules.)

### C. Channel simplification

- [x] **C1. Replace the nine re-key → combine → strip patterns with `pairById`** or with id-keyed
  channels: `preprocessing.nf` (6), `projection.nf` (2; updated by D27d-3), `benchmarking.nf` (1),
  `curation.nf` (1). Build `FASTA_SIZES` pairing once.
- [x] **C2. One gating mechanism for optional branches (P10).** (updated by D27e-2) (updated by D27a-4: the V1 log
  shows the warning for the `if`-gated decoy processes with fully qualified selectors. Re-check
  it after 27d's `.*:NAME$` selectors.)
  - Test whether a `withName` selector whose process is included but never invoked (because of
    an `if`) raises Nextflow's "no process matches selector" warning.
  - If it does not, gate curation with `if (params.curate)` in `fomo.nf`, like the decoy branch,
    and drop `ext.when`.
  - If it does, keep `ext.when` for curation and record why.
  - Remove the `when:` block from CURATE_MODELS if it becomes unused.
- [x] ~~**C3. Move the MultiQC input tuple builder, duplicated in `reporting.nf` and
  `run_summary.nf`, into a helper in `fomo_utils.nf`**~~ (or delete it if 27i's design makes it
  moot; note it). → Moot: since 27c each builder is one `tuple(meta, files, [two configs])`, and
  27i B3 redesigns the MULTIQC inputs.
- [x] **C4. Nextflow lint warnings → 0** (unused params, implicit `it`, deprecated `Channel`).
  Tighten the pre-commit hook so warnings fail too.

## Acceptance

- G1, G2 (lint with 0 warnings), G3 (stub task counts unchanged), and G4 identical.
- The A3 grep returns nothing. `meta.gtype` does not exist.
- The V2 trace has unique tags per process (B2).
- Results record the line counts of the subworkflows and `modules.config`, before and after.

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|
| D27e-1 | A1 | helpers `trackOf`, `withTrack`, `pairById`, `asList` | `trackOf`, `withTrack`, `pairBy(left, right, key)` + `pairById` (key `meta.id`); no `asList` | the strict syntax allows neither default parameters (`key = {…}`) nor calling a closure variable as `key(x)` (`nextflow lint`: "`key` is not defined"), hence two functions and `key.call(…)`. No `instanceof List` site is left after 27b; the remaining `[x].flatten()` idiom is already one expression | none |
| D27e-2 | C2 | if `if`-gating warns, keep `ext.when` for curation | both optional branches are gated by emptying their input channel (`params.include_decoy ? ch_sources : channel.empty()`, `params.curate ? ch_allmodels_raw : channel.empty()`); the processes are always invoked and run no task when off. No `when:` in any module, no `ext.when` | `if`-gating does warn: the 27d V1 log has "There's no process matching config selector: .*:BEDTOOLS_COMPLEMENT$" (and EXTRACT_DECOY_SEQUENCES, RELOCATE_LOCI). Keeping `ext.when` for curation would leave two mechanisms and the decoy warnings; an empty input gives one mechanism, no warnings (the S1 stub log has no WARN line) and unchanged task counts (G3) | master §8 "`when:` only where config gates the module" → no module has `when:` (enforced by `tests/test_repo_conventions.py`); 27j documents the gating rule |

## Results

- **Track keys** are stamped by `withTrack` where a feature type or the decoy flag is first set
  (the PREPROCESSING candidates, RELOCATE_LOCI's input) and rebuilt from them in PROJECTION's
  merge. `grep -rn "meta.decoy ? '" workflows subworkflows modules conf` and `grep -rn gtype …`
  return nothing; TD2 branches on `meta.td2`; RENAME_FASTA_HEADERS writes `meta.label`; the
  sample-name ladders are `"${order}_${label}"` (+1 after `1_raw` on the annotation-stats rows).
  The renamed intermediate is `allsources.<track>.fasta` (decoys: `.decoy` instead of `_decoy`).
- **Identity:** target-scoped tasks carry `id = target_id` (no `'allModels'` or `'curated'` meta
  ids); every module but RUN_SUMMARY_TABLES uses the master §8 tag rule (enforced by
  `tests/test_repo_conventions.py`). In the V2 trace no two tasks of a process share a tag;
  MINIMAP2_ALIGN shows `<target>.<track>`.
- **Pairing:** 9 re-key → combine → strip patterns became `pairById`/`pairBy` calls
  (`preprocessing.nf` 5, `projection.nf` 2, `benchmarking.nf` 1, `curation.nf` 1); the FASTA
  sizes are keyed once (`ch_sizes`). The FILTER_ANNOTATION input stays a `join(…,
  failOnMismatch: true)`: it must fail, not drop, on a missing partner. Process outputs are
  immutable tuples, so `pairBy` builds a new list (`item.toList() + [value]`; `item + [value]`
  threw `UnsupportedOperationException` in the first G3 attempt).
- **Gating (D27e-2):** decoys and curation are gated by an empty input channel; the V1 and V2
  logs have no warning at all. EXTRACT_SEQUENCES and EXTRACT_DECOY_SEQUENCES share one selector
  (`.*:EXTRACT_(DECOY_)?SEQUENCES$`, prefix `<id>.<track>.spliced`).
- **Lint:** 0 errors, 0 warnings; the pre-commit hook runs `dev/nextflow_lint.sh`, which fails on
  warnings (checked with a throw-away file that declares an unused variable).
- **Line counts** (before → after): `preprocessing.nf` 172 → 131, `projection.nf` 113 → 106,
  `curation.nf` 29 → 28, `benchmarking.nf` 57 → 57, `fomo_utils.nf` 19 → 52 (the helpers),
  `conf/modules.config` 272 → 222, modules 964 → 952.
- **Gates:** G1 122 passed (incl. the run-summary fixture test, so the sample names are
  byte-identical); G2 pass; G3 S1–S4 pass, task counts unchanged; G4 V1–V4 identical. Cluster use
  0.4 CPU-h.
