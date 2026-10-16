# 27j — Documentation, plan records, dev tools, comment budget

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Review items: V5, V6, V7, T4,
T5, D5, and F16/F17/F19/F23 of the docs review.

## Status: done (2026-10-15)

## Goal

Each fact has **one home**:

- **code and tests:** behaviour and invariants;
- **one short comment:** why;
- **CLAUDE.md:** an index of the invariants and conventions, with checklists;
- **README and `docs/`:** for people;
- **`plans/` with a status:** decisions;
- **CHANGELOG:** history.

The comment budget is enforced mechanically from now on.

## Steps

### A. CLAUDE.md (≤ 220 lines)

- [x] **A1. Rewrite CLAUDE.md from the code as it is now, not by trimming the old file.**
  Sections:
  1. Project context: what FOMO does, inputs and outputs.
  2. **Development guidance, updated for U2.** Use the Seqera MCP tools to find nf-core prior
     art, then vendor the module into `modules/local`, trimmed to what FOMO calls. The meta-map
     and `ext.args` conventions stay. Wave containers are declared in `conf/containers.config`.
     Validate config before launching.
  3. Glossary (target, source, both, track, record) and roles.
  4. **Invariants.** Each is one or two present-tense lines plus the test or code that
     enforces it (master §4).
  5. Architecture map: subworkflows, the meta keys, the data records, and the publish layout
     (pointing to `docs/output.md`).
  6. Conventions: comments, naming, CLI, containers, config (master §8).
  7. **Checklists:** add a track; add a stats section; add a param; add a module.
  8. How to test: pytest, nf-test, golden, CI.
- [x] **A2. Leave out:**
  - history, incidents, measurements and run names;
  - the legacy AGAT table;
  - per-plan narratives;
  - anything the code or tests already state.

### B. People-facing documentation

- [x] **B1. README** as the entry point:
  - what FOMO is;
  - requirements (Nextflow version floor, Singularity or Docker);
  - quick start (`-profile test,singularity` or `test,docker`);
  - a samplesheet example and the role semantics (once);
  - the key params (pointing to `--help`);
  - outputs (pointing to `docs/output.md`);
  - how to cite and the licence.
- [x] **B2. `docs/usage.md`, `docs/output.md` (the publish tree, once),
  `docs/architecture.md`, and `docs/data_records.md`** (the `fomo-data/1` schema per section,
  generated from or checked against `fomo_stats.SECTIONS`).

### C. Plans as decision records

- [x] **C1. Give every existing plan (00–26) an ADR header.** It has `Status:` (Implemented,
  Superseded by NN, or Partially superseded by NN), a `Date`, and a 3-line decision summary.
  **Only the header is added; the bodies stay untouched.**
- [x] **C2. Write `plans/README.md`** (updated by D27j-1), an index with one line per plan and its status. Mark
  27–27k and the frozen `27_refactor_original/`.

### D. History, contributing and dev tools

- [x] **D1. `CHANGELOG.md`.**
  - `0.2.0` (this refactor) lists the user-visible changes: the species pattern, timestamped
    `pipeline_info`, MultiQC changes (software versions, no duplicate target rows, renderer),
    the conda profile removed, modules local, pinned images.
  - Summarise the earlier history in one `0.1.0` entry, which includes the "results not
    comparable before plans/20" note.
  - Bump `manifest.version` to `0.2.0dev`.
- [x] **D2. `CONTRIBUTING.md`:** the principles (review §4), conventions, gates, the commit
  format, and how to run CI locally.
- [x] **D3. Create `dev/` for maintenance tools** (updated by D27j-4) (it already holds `environment.yml`, 27a, and
  `nextflow_lint.sh`, 27e C4). Move `bin/subsample_test_data.sh` to `dev/`
  (fix its header and paths). Move `plans/13_resource_rightsizing_audit.py` to `dev/` and leave
  a pointer in plans/13. Copy `ablation/scripts/code_metrics.py` to `dev/code_metrics.py`; the
  ablation copy stays as the study record.
- [x] **D4. Budget mode.** (updated by D27j-2) (updated by D27b-4: match the narrower phrase lists of 27b's C3/C4
  greps, not the review's regexes, which flag present-tense phrases such as "was dropped" and
  "verified against".) Add `dev/code_metrics.py --check` with thresholds:
  - Nextflow layer comment ratio ≤ 0.40;
  - 0 incident hits; 0 history hits;
  - no `plans/` reference in user-facing files.

  Run it in pre-commit and CI.

### E. Checks

- [x] **E1. Link check:** every relative link and path mentioned in CLAUDE.md, README, `docs/`
  and `plans/README.md` exists. Write a small script in `dev/` and run it in CI.
- [x] **E2. Duplication check:** each of the D5 topics (role semantics, publish tree, unique IDs,
  gzip/faidx, trans-spliced genes) is explained in exactly one place. Grep for its key phrase
  and record the result in Results.

## Acceptance

- G1, G2 (now including the budget and the link check), G3, and G4 identical. Docs do not change
  outputs.
- CLAUDE.md is ≤ 220 lines and every invariant line names its enforcing test.
- Results record the master §2 metric table, before and after.

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|
| D27j-1 | C2 | write `plans/README.md` | generated by `dev/plans_index.py` from each plan's header or `## Status:` line | a hand-written index goes stale with every status change; plan titles are written without backticks (one names a script that no longer exists) | 27k reruns `dev/plans_index.py` |
| D27j-2 | D4 | no `plans/` reference in user-facing files | checked in README, `docs/`, both schemas and `assets/multiqc/`; report text written by `bin/` scripts is not parsed. The one such reference (`(plans/23)` in the run-summary heatmap description) was removed by hand, and the six heatmap plot specs of the run-summary fixtures changed in that description only | parsing strings out of Python for a one-off is not worth it | none |
| D27j-3 | A1 (invariant 8) | each invariant names its enforcing test (master §4 #8: "27a pytest") | `tests/modules/test_td2_predict.py` added: empty FASTA, no complete ORF, normal case | there was no test of the TD2 guards; 27a did not write the one §4 assigned to it | none |
| D27j-4 | D3 | move `plans/13_resource_rightsizing_audit.py` to `dev/` | `dev/resource_audit.py` (report printing split into two functions for the statement limit); `dev/subsample_test_data.sh` takes the downloads directory as an argument | `dev/` is linted, `plans/` is not; the old script assumed an `annotation_downloads/` directory inside the repository, which does not exist | none |

## Results

- **CLAUDE.md** rewritten from the code: 514 → 150 lines; 16 invariants, each naming what
  enforces it; checklists for a track, a statistics section, a parameter and a module.
- **People-facing docs:** README (entry point, the samplesheet and roles once), `docs/usage.md`,
  `docs/output.md` (the publish tree once), `docs/architecture.md`, `docs/data_records.md`
  (checked against `fomo_stats.SECTIONS` by `tests/test_docs.py`); CHANGELOG (0.2.0, 0.1.0) and
  CONTRIBUTING; `manifest.version` 0.2.0dev.
- **Plans:** decision-record headers on the 30 plans 00–26 (bodies untouched), `plans/README.md`
  generated.
- **dev/:** `code_metrics.py --check` (Nextflow workflows 0.26, config 0.20, budget 0.40; no
  history or run phrase; no `plans/` in user-facing files), `check_links.py`, `plans_index.py`,
  `resource_audit.py`, `subsample_test_data.sh`, tested in `tests/dev`; the budget and link checks
  run in pre-commit and so in CI.
- **E2, one home per topic** (grep over CLAUDE.md, README, CONTRIBUTING, CHANGELOG, `docs/`):
  role semantics — README only; publish tree — `docs/output.md` only (`docs/usage.md` points at
  two paths in it); unique projected IDs — explained in `docs/architecture.md`, stated as a file
  format in `docs/output.md` and as invariant 2 in CLAUDE.md; gzip/faidx — explained in
  `docs/architecture.md`, invariant 4 states the rule; trans-spliced models — explained only in
  `bin/filter_transcript.py`'s docstring (and the MultiQC column description).
- **Master §2 metrics** (before → now):

  | metric | before | now |
  |:---|---:|---:|
  | comment : code, Nextflow workflows / config | 0.83 / 0.88 | 0.26 / 0.20 |
  | run / incident / history phrases (budget lists) | 78 (review regexes) | 0 (review regexes: 5, all present-tense "measured"/"verified" or the variable `observed`) |
  | `CLAUDE.md` | 514 lines | 150 lines |
  | copies of the track-name expression | 35 | 0 |
  | parsers of the per-source filename grammar | 3 + 1 partial | 0 |
  | MultiQC adapter processes | 5 | 0 |
  | GFF statistics modules | 3 | 1 |
  | `modules/nf-core` | 9 modules, 2 forked | gone |
  | tests | 1 script | 275 pytest, 4 nf-test scenarios, golden V1–V4, CI |
  | `nextflow lint` errors / warnings | 6 / — | 0 / 0 |
  | Python functions above complexity 15 | 6 | 0 (maximum 12) |
  | tasks in V2 | 167 | 133 |
- **Gates:** G1 275 passed; G2 pass (now with the comment budget and the link check); G3 S1–S4
  pass; G4 V1–V4 identical. Cluster use 0.3 CPU-h.
