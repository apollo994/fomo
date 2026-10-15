# 27d — Config structure: include order, selectors, resources, site profile, pipeline_info

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decisions: U5 (timestamped
`pipeline_info`), T5. Review items: P5, P6, C4, R7, F13/F14 of the workflow review, N5 (tabs).

## Status: done (2026-10-15)

## Goal

Configuration that can be read top to bottom with no precedence traps:

- `base.config` holds every resource, with labels that state real needs;
- `crg.config` holds only site settings;
- `containers.config` holds every image (from 27c);
- `modules.config` holds only `ext.*` and `publishDir`, with one selector style.

Process placement in the call tree no longer matters.

## Steps

- [x] **A1. Record the resolved config before changing anything.**
  `nextflow config -flat -profile crg`, `-profile test,singularity` and `-profile docker`
  → `/nfs/scratch01/rg/fzanarello/fomo_refactor/config_before/`.
- [x] **A2. Include order.** In `nextflow.config`:
  1. params
  2. `includeConfig 'conf/base.config'`
  3. `includeConfig 'conf/containers.config'`
  4. `profiles { … }`
  5. `includeConfig 'conf/modules.config'`

  Delete the "later selector wins / precedence trap" comments (`base.config`, `crg.config`).
- [x] **A3. Site profile.** `crg.config` holds settings only, with no `profiles { crg { … } }`
  wrapper; `nextflow.config` wraps it the way it wraps `test.config`. It keeps:
  - executor, queue, `clusterOptions` (account), `cache = 'lenient'`, `resourceLimits`;
  - executor tuning (`pollInterval`, `queueSize`, `dumpInterval`);
  - `singularity.*`;
  - the MINIMAP2_INDEX `maxForks` (a site QOS limit).

  `workDir = System.getenv('NXF_WORK') ?: '/nfs/scratch01/rg/fzanarello/work_fomo'`, so the
  default stays as it is. Replace tabs with spaces in every config file.
- [x] **A4. Resources live only in `base.config`.** Labels that state real needs:

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
- [x] **A5. Selectors.** (updated by D27d-3)
  - Every `withName` in `modules.config` becomes `'.*:NAME$'`.
  - The two TD2_NONCODING instances are told apart by `meta.stage` (`'input'` or
    `'projected'`), set by the callers in `preprocessing.nf` and `projection.nf`.
  - The ext closures of `TD2_CODING_FILTER` read `meta.stage`. TD2_PREDICT keeps one block,
    because both instances share the same args.
  - Check: `nextflow config` shows each process resolving the same `ext.*` as before (diff of
    A1's flat dump against the new one, keys `process.$NAME.ext*`).
- [x] **A6. Timestamped `pipeline_info` (U5).** (updated by D27d-1) `trace`, `report`, `timeline` and `dag` files get a
  run timestamp in the name: a hidden param `trace_report_suffix`, default
  `new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')` (D27d-1). Drop `overwrite = true`. Fix the
  comment.
- [x] **A7. Manifest and schema.** (updated by D27d-2)
  - `homePage 'https://github.com/apollo994/fomo'`, `defaultBranch 'main'`, an author (as
    `contributors`, D27d-2).
  - `nextflowVersion`: raise to `!>=25.04.0` (updated by D27c-6: necessary — topic channels
    leave preview in 25.04).
  - Fix the `$id` of `nextflow_schema.json` and `assets/schema_input.json`
    (`https://raw.githubusercontent.com/apollo994/fomo/main/...`).
- [x] **A8. Delete the comments these changes made obsolete** (master §8).
- [x] **A9. Compare the resolved config again (A1).** Every difference must be intended. List them
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
| D27d-1 | A6 | `def trace_timestamp = …` in `nextflow.config` | hidden param `params.trace_report_suffix` (default: the launch time), added to `nextflow_schema.json` | the strict config syntax rejects variable declarations next to config statements (`nextflow lint`: "Variable declarations cannot be mixed with config statements", hit in 27c); a param is the form nf-core uses and lets a user name a run's reports | 27j: `docs/usage.md` may mention the hidden param |
| D27d-2 | A7 | an `author` in the manifest | `manifest.contributors` (name, affiliation, GitHub URL, roles; no e-mail) | `manifest.author` is deprecated since Nextflow 24.10 in favour of `contributors` (Nextflow manifest reference, docs.seqera.io/nextflow/reference/config/manifest). `defaultBranch` is kept as planned; that reference deprecates it only from 26.04 | none |
| D27d-3 | A5 | the TD2_CODING_FILTER ext closures read `meta.stage` | they do (prefix, sample name); the module reads `meta.stage` itself for `--stage`, so `ext.stage` is gone. In `projection.nf` the coding-id join is keyed on `target_id` | the TD2 output meta now carries `stage`, the BAM_TO_GFF meta does not, so a whole-meta join would no longer match (one lncRNA drop list per target, so `target_id` is the key) | 27e C1: `projection.nf` has two keyed pairings, not one (updated by D27d-3) |
| D27d-4 | A1, A9 | compare `nextflow config -flat` dumps | the flat dumps are kept (and the `ext.*` keys diffed), and the per-process comparison uses the effective values: a stub run under `-profile test,crg` with the local executor and the V2 flags (all 31 processes), cpus/memory/time/container from its trace (`fomo_refactor/resolve.sh`) | the flat dump lists the declarations, not what each process resolves to; the stub run resolves them the way the real run does | none |

## Results

**A9, `-profile crg` (the production profile):** for all 31 processes (V2 flags) the effective
cpus, memory, time and container are identical before and after
(`fomo_refactor/config_{before,after}/stub_crg/resolved.tsv`). MINIMAP2_INDEX (3 CPUs, 32 GB × attempt,
1 h × attempt, `maxForks 96`), MINIMAP2_ALIGN (2 CPUs, 24 GB × 2^(attempt−1), 1 h × attempt,
`maxRetries 2`) and RUN_SUMMARY_TABLES (8 GB × attempt, 2 h × attempt) are unchanged; the two batch
processes keep 4 CPUs, 4 GB × attempt, 2 h × attempt and the 123/125 retry codes. `workDir` is
unchanged unless `NXF_WORK` is set.

**Intended differences outside `crg`** (A4 moves the genome-size-driven values into
`base.config`):

| process | before (`test`, `docker`) | after |
|:---|:---|:---|
| MINIMAP2_INDEX | `process_low`: 1 CPU, 4 GB × attempt, 15 min × attempt | `process_index`: 3 CPUs, 32 GB × attempt, 1 h × attempt (`test` caps: 2 CPUs, 4 GB, 30 min) |
| MINIMAP2_ALIGN | 2 CPUs, 24 GB × attempt, 1 h × attempt, `maxRetries 3` | `process_align`: 24 GB × 2^(attempt−1), `maxRetries 2` |
| RUN_SUMMARY_TABLES | `process_single`: 2 GB × attempt, 10 min × attempt | 8 GB × attempt, 2 h × attempt (`test` caps: 4 GB, 30 min) |

**`ext.*`:** the only differences in the three flat dumps are the TD2_CODING_FILTER blocks (two
merged into one keyed on `meta.stage`, `ext.stage` removed). Every `withName` in
`modules.config` is `.*:NAME$` (38 selectors). `process_medium` and `process_high` are gone.

**Other:** `nextflowVersion '!>=25.04.0'` (D27c-6), `homePage`, `defaultBranch`, `contributors`;
both schema `$id`s point at `apollo994/fomo/main`; no tabs in any config file; `--help` renders.

**pipeline_info:** two consecutive runs into one outdir keep both sets
(`execution_{report,timeline,trace}_<time>`, `pipeline_dag_<time>`).

**Gates:** G1 122 passed; G2 pass (`nextflow lint` 43 files, no errors or warnings); G3 S1–S4
pass; G4 V1–V4 identical (`gate/27d`; the V1 `pipeline_info/` holds timestamped files).
Cluster use 0.5 CPU-h.
