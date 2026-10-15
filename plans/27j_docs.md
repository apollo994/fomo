# 27j — Documentation, plan records, dev tools, comment budget

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Review items: V5, V6, V7, T4,
T5, D5, and F16/F17/F19/F23 of the docs review.

## Status: not started

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

- [ ] **A1. Rewrite CLAUDE.md from the code as it is now, not by trimming the old file.**
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
- [ ] **A2. Leave out:**
  - history, incidents, measurements and run names;
  - the legacy AGAT table;
  - per-plan narratives;
  - anything the code or tests already state.

### B. People-facing documentation

- [ ] **B1. README** as the entry point:
  - what FOMO is;
  - requirements (Nextflow version floor, Singularity or Docker);
  - quick start (`-profile test,singularity` or `test,docker`);
  - a samplesheet example and the role semantics (once);
  - the key params (pointing to `--help`);
  - outputs (pointing to `docs/output.md`);
  - how to cite and the licence.
- [ ] **B2. `docs/usage.md`, `docs/output.md` (the publish tree, once),
  `docs/architecture.md`, and `docs/data_records.md`** (the `fomo-data/1` schema per section,
  generated from or checked against `fomo_stats.SECTIONS`).

### C. Plans as decision records

- [ ] **C1. Give every existing plan (00–26) an ADR header.** It has `Status:` (Implemented,
  Superseded by NN, or Partially superseded by NN), a `Date`, and a 3-line decision summary.
  **Only the header is added; the bodies stay untouched.**
- [ ] **C2. Write `plans/README.md`**, an index with one line per plan and its status. Mark
  27–27k and the frozen `27_refactor_original/`.

### D. History, contributing and dev tools

- [ ] **D1. `CHANGELOG.md`.**
  - `0.2.0` (this refactor) lists the user-visible changes: the species pattern, timestamped
    `pipeline_info`, MultiQC changes (software versions, no duplicate target rows, renderer),
    the conda profile removed, modules local, pinned images.
  - Summarise the earlier history in one `0.1.0` entry, which includes the "results not
    comparable before plans/20" note.
  - Bump `manifest.version` to `0.2.0dev`.
- [ ] **D2. `CONTRIBUTING.md`:** the principles (review §4), conventions, gates, the commit
  format, and how to run CI locally.
- [ ] **D3. Create `dev/` for maintenance tools** (it already holds `environment.yml`, 27a, and
  `nextflow_lint.sh`, 27e C4). Move `bin/subsample_test_data.sh` to `dev/`
  (fix its header and paths). Move `plans/13_resource_rightsizing_audit.py` to `dev/` and leave
  a pointer in plans/13. Copy `ablation/scripts/code_metrics.py` to `dev/code_metrics.py`; the
  ablation copy stays as the study record.
- [ ] **D4. Budget mode.** (updated by D27b-4: match the narrower phrase lists of 27b's C3/C4
  greps, not the review's regexes, which flag present-tense phrases such as "was dropped" and
  "verified against".) Add `dev/code_metrics.py --check` with thresholds:
  - Nextflow layer comment ratio ≤ 0.40;
  - 0 incident hits; 0 history hits;
  - no `plans/` reference in user-facing files.

  Run it in pre-commit and CI.

### E. Checks

- [ ] **E1. Link check:** every relative link and path mentioned in CLAUDE.md, README, `docs/`
  and `plans/README.md` exists. Write a small script in `dev/` and run it in CI.
- [ ] **E2. Duplication check:** each of the D5 topics (role semantics, publish tree, unique IDs,
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

## Results

*(fill in)*
