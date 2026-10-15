# 27d — Config structure: include order, selectors, resources, site profile, pipeline_info

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decisions: U5 (timestamped
`pipeline_info`), T5. Review items: P5, P6, C4, R7, F13/F14 of the workflow review, N5 (tabs).

## Status: not started

## Goal

Configuration that can be read top to bottom with no precedence traps:

- `base.config` holds every resource, with labels that state real needs;
- `crg.config` holds only site settings;
- `containers.config` holds every image (from 27c);
- `modules.config` holds only `ext.*` and `publishDir`, with one selector style.

Process placement in the call tree no longer matters.

## Steps

- [ ] **A1. Record the resolved config before changing anything.**
  `nextflow config -flat -profile crg`, `-profile test,singularity` and `-profile docker`
  → `/nfs/scratch01/rg/fzanarello/fomo_refactor/config_before/`.
- [ ] **A2. Include order.** In `nextflow.config`:
  1. params
  2. `includeConfig 'conf/base.config'`
  3. `includeConfig 'conf/containers.config'`
  4. `profiles { … }`
  5. `includeConfig 'conf/modules.config'`

  Delete the "later selector wins / precedence trap" comments (`base.config`, `crg.config`).
- [ ] **A3. Site profile.** `crg.config` holds settings only, with no `profiles { crg { … } }`
  wrapper; `nextflow.config` wraps it the way it wraps `test.config`. It keeps:
  - executor, queue, `clusterOptions` (account), `cache = 'lenient'`, `resourceLimits`;
  - executor tuning (`pollInterval`, `queueSize`, `dumpInterval`);
  - `singularity.*`;
  - the MINIMAP2_INDEX `maxForks` (a site QOS limit).

  `workDir = System.getenv('NXF_WORK') ?: '/nfs/scratch01/rg/fzanarello/work_fomo'`, so the
  default stays as it is. Replace tabs with spaces in every config file.
- [ ] **A4. Resources live only in `base.config`.** Labels that state real needs:

  | label | cpus | memory | time | notes |
  |:---|---:|:---|:---|:---|
  | `process_single` | 1 | 2 GB × attempt | 10 min × attempt | |
  | `process_low` | 1 | 4 GB × attempt | 15 min × attempt | |
  | `process_batch` | 4 | 4 GB × attempt | 2 h × attempt | the two batch modules; replaces `process_high` plus its override |
  | `process_index` | 3 | 32 GB × attempt | 1 h × attempt | MINIMAP2_INDEX; moved from `crg.config`, it is genome-size driven |
  | `process_align` | 2 | 24 GB × 2^(attempt−1) | 1 h × attempt | MINIMAP2_ALIGN, `maxRetries = 2`; moved from `crg.config` + `base.config` |

  - RUN_SUMMARY_TABLES (8 GB × attempt, 2 h × attempt) moves from `crg.config` into a
    `withName` block in `base.config`.
  - Delete `process_medium` and `process_high` once no module uses them.
  - Keep the retry list and the 123/125 additions from 27b.
- [ ] **A5. Selectors.**
  - Every `withName` in `modules.config` becomes `'.*:NAME$'`.
  - The two TD2_NONCODING instances are told apart by `meta.stage` (`'input'` or
    `'projected'`), set by the callers in `preprocessing.nf` and `projection.nf`.
  - The ext closures of `TD2_CODING_FILTER` read `meta.stage`. TD2_PREDICT keeps one block,
    because both instances share the same args.
  - Check: `nextflow config` shows each process resolving the same `ext.*` as before (diff of
    A1's flat dump against the new one, keys `process.$NAME.ext*`).
- [ ] **A6. Timestamped `pipeline_info` (U5).** `trace`, `report`, `timeline` and `dag` files get a
  run timestamp in the name. Use one `def trace_timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')`
  in `nextflow.config`, and drop `overwrite = true`. Fix the comment.
- [ ] **A7. Manifest and schema.**
  - `homePage 'https://github.com/apollo994/fomo'`, `defaultBranch 'main'`, an `author`.
  - `nextflowVersion`: raise to `!>=25.04.0` if 27c-D4 found it necessary.
  - Fix the `$id` of `nextflow_schema.json` and `assets/schema_input.json`
    (`https://raw.githubusercontent.com/apollo994/fomo/main/...`).
- [ ] **A8. Delete the comments these changes made obsolete** (master §8).
- [ ] **A9. Compare the resolved config again (A1).** Every difference must be intended. List them
  in Results: process → cpus/memory/time/maxRetries/container, before and after, for `crg`. The
  `crg` values for MINIMAP2_INDEX, MINIMAP2_ALIGN and RUN_SUMMARY_TABLES must be unchanged.

## Acceptance

- G1, G2 (`nextflow lint` clean, with no new warnings), G3, and G4 identical.
- The A9 diff is fully explained in Results.
- `pipeline_info/` files of a G4 run are timestamped. Two consecutive runs into the same outdir
  keep both traces.

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|

## Results

*(fill in)*
