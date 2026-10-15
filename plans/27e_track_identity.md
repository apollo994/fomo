# 27e — Track registry, meta identity, channel simplification

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decision: T4. Review items:
N1, N2, P7, P10, and F3/F5/F10/F11/F23/F25 of the workflow review.

## Status: not started

## Goal

A track — feature type plus decoy — is defined **once** and carried in meta. No expression
recomputes it, no task is tagged `allModels`, and channel pairing uses one helper instead of nine
hand-written key/strip patterns.

## Steps

### A. Registry and meta keys

- [ ] **A1. Write `subworkflows/local/fomo_utils.nf`** (it may already exist from 27c-D2) with:
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
- [ ] **A2. Stamp the track keys once.**
  - PREPROCESSING: where `feature_type`/`decoy` are first set, for the real track
    (`preprocessing.nf` ≈ line 50) and for the decoy track (≈ line 126).
  - PROJECTION's merge: replace `gtype` with the carried keys.
  - Delete `meta.gtype` (underscore form). MERGE_SOURCE_FASTA's prefix becomes
    `allsources.${meta.track}`; that file is intermediate, so U1 is not affected.
- [ ] **A3. Replace every `${meta.feature_type}${meta.decoy ? '.decoy' : ''}` with
  `${meta.track}`** (35 sites: `modules.config`, modules, subworkflows). Check with
  `grep -rn "meta.decoy ? '" workflows subworkflows modules conf`, which must return nothing.
- [ ] **A4. TD2 branching uses `meta.td2`** (`preprocessing.nf`, `projection.nf`), not
  `!meta.decoy && meta.feature_type == 'lnc_RNA'`.
- [ ] **A5. The sample-name ladders in `modules.config` use `meta.order`/`meta.label`** for
  tracks. The kind ladder for raw/filtered stays until 27i deletes these closures. **Sample
  names must stay byte-identical**, because `run_summary_tables.py` parses them until 27h:
  G4 + the run-summary fixture test.

### B. `meta.id` and tags

- [ ] **B1. Stop setting `meta.id` to pseudo values (`'allModels'`, `'curated'`).**
  - Target-scoped tasks carry `id = target_id`.
  - Every prefix that relied on the pseudo id is spelled out with `meta.target_id`/`meta.track`
    in `modules.config` (most already are).
  - The pseudo-source `allModels_raw` in the stats grammar is untouched; 27g replaces it.
- [ ] **B2. One tag rule** (master §8) in every module. Check: in the V2 trace, no two tasks of
  the same process share a tag. MINIMAP2_ALIGN shows `<target>.<track>`.

### C. Channel simplification

- [ ] **C1. Replace the nine re-key → combine → strip patterns with `pairById`** or with id-keyed
  channels: `preprocessing.nf` (6), `projection.nf` (1), `benchmarking.nf` (1),
  `curation.nf` (1). Build `FASTA_SIZES` pairing once.
- [ ] **C2. One gating mechanism for optional branches (P10).**
  - Test whether a `withName` selector whose process is included but never invoked (because of
    an `if`) raises Nextflow's "no process matches selector" warning.
  - If it does not, gate curation with `if (params.curate)` in `fomo.nf`, like the decoy branch,
    and drop `ext.when`.
  - If it does, keep `ext.when` for curation and record why.
  - Remove the `when:` block from CURATE_MODELS if it becomes unused.
- [ ] **C3. Move the MultiQC input tuple builder, duplicated in `reporting.nf` and
  `run_summary.nf`, into a helper in `fomo_utils.nf`** (or delete it if 27i's design makes it
  moot; note it).
- [ ] **C4. Nextflow lint warnings → 0** (unused params, implicit `it`, deprecated `Channel`).
  Tighten the pre-commit hook so warnings fail too.

## Acceptance

- G1, G2 (lint with 0 warnings), G3 (stub task counts unchanged), and G4 identical.
- The A3 grep returns nothing. `meta.gtype` does not exist.
- The V2 trace has unique tags per process (B2).
- Results record the line counts of the subworkflows and `modules.config`, before and after.

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|

## Results

*(fill in)*
