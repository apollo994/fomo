# 27k — Final validation, real-data regression, plan comparison

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first).

## Status: not started

## Goal

Prove the refactor met the plan:

- every gate passes on a clean checkout;
- the real-data deliverables match the pre-refactor run;
- the goals in master §2 are measured.

Then **compare the working plan with the frozen original** and validate every deviation.

## Steps

### A. Clean-checkout gates

- [ ] **A1. Set up a clean worktree.** `git worktree add` a clean checkout of `refactor`
  outside the repo (scratch root), then create the dev env there from `dev/environment.yml`.
- [ ] **A2. Run G1, G2, G3 and G4 there, with every variant run fresh** (no `-resume`). Run every
  CI `run:` command locally with Singularity in place of Docker. List what could not be
  verified without pushing.

### B. Real-data regression (G5)

- [ ] **B1. Run the 23 Lycaenidae `both` species with `--include_decoy`, `-profile crg`.**
  - Samplesheet: `tests/golden/README.md` (from 27a-A8).
  - Use a dedicated `-w /nfs/scratch01/rg/fzanarello/fomo_refactor/work_g5` and outdir
    `…/fomo_refactor/g5/`.
  - Launch the head job as in the ablation study: `sbatch`, QOS `pipelines`, 1 cpu, 4 GB,
    12 h limit.
- [ ] **B2. Snapshot the outputs and `compare.py` them with `tests/golden/realdata_B0.tsv`.**
  **Zero differences** is required. A difference is a stop condition (master §7c): investigate
  it, never accept it silently.
- [ ] **B3. Record CPU-h used** in master §10. The total must stay within the budget of stop
  condition (d).

### C. Measure the goals

- [ ] **C1. Fill the master §2 table with "after" values:**
  - `dev/code_metrics.py`;
  - grep counts;
  - `ruff --statistics`;
  - `nextflow lint`;
  - task counts from the V2 trace (27a-A9 against now);
  - the line counts per layer;
  - the length of CLAUDE.md.

  Every miss gets a sentence.

### D. Plan comparison and deviation validation

- [ ] **D1. Diff the original against the working plan:** `diff -ru plans/27_refactor_original/
  plans/ --include` the `27*.md` files only. Save the raw diff to
  `/nfs/scratch01/rg/fzanarello/fomo_refactor/plan_diff.txt`.
- [ ] **D2. Classify every hunk:**
  - **bookkeeping:** status, checkboxes, Results, the status log;
  - **documented deviation:** it has a `D27x-N` entry;
  - **undocumented change.**

  Every undocumented change gets a retroactive `D27x-N` entry, or is reverted if it was a
  mistake.
- [ ] **D3. Validate each deviation in `plans/27z_deviation_review.md`**, in a table:

  ```
  | id | subplan | summary | class (correction / improvement / scope cut / scope addition) | consistent with U1–U8? | invariants and gate preserved? | evidence present? | verdict |
  ```

  - A deviation is **valid** only if it keeps every U-decision, keeps the gate and the
    invariants, and its evidence is reproducible (command, output or test).
  - **Invalid deviations are fixed:** change the code back toward the plan or complete the
    missing work, re-run the gates, and record it.
  - A deviation that needs a U-decision change has already been escalated (master §7). Confirm
    the user's answer is recorded.
- [ ] **D4. Plan coverage:** every step of every original subplan is done, or appears in a
  validated deviation as replaced or cut. List any that were silently skipped (there must be
  none).

### E. Close-out

- [ ] **E1. Write master §11 "Outcome":**
  - the goals table (before and after);
  - G5 result;
  - task-count change;
  - lines removed;
  - the deviation summary (counts by class, all valid);
  - what CI could not verify;
  - follow-ups. These include the scientific recommendations from
    `ablation/ABLATION_REPORT.md` (adaptive donor set, per-track thresholds, mRNA coding mask,
    `splice:hq` option, decoy FDR estimator) as candidate plans 28+.
- [ ] **E2. Final commit:** `refactor(27k): final validation and plan comparison`. Set the master
  status to `done (YYYY-MM-DD)`. Do not push; tell the user the branch is ready for review.

## Acceptance

- Every gate passes on the clean checkout.
- G5 has zero differences.
- `27z_deviation_review.md` has no invalid deviation.
- Master §11 is written.

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|

## Results

*(fill in)*
