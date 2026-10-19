# 27k — Final validation, real-data regression, plan comparison

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first).

## Status: done (2026-10-16)

## Goal

Prove the refactor met the plan:

- every gate passes on a clean checkout;
- the real-data deliverables match the pre-refactor run;
- the goals in master §2 are measured.

Then **compare the working plan with the frozen original** and validate every deviation.

## Steps

### A. Clean-checkout gates

- [x] **A1. Set up a clean worktree.** (updated by D27k-2) `git worktree add` a clean checkout of `refactor`
  outside the repo (scratch root), then create the dev env there from `dev/environment.yml`.
- [x] **A2. Run G1, G2, G3 and G4 there, with every variant run fresh** (no `-resume`). Run every
  CI `run:` command locally with Singularity in place of Docker. List what could not be
  verified without pushing.

### B. Real-data regression (G5)

- [x] **B1. Run the 23 Lycaenidae `both` species with `--include_decoy`, `-profile crg`.**
  - Samplesheet: `tests/golden/README.md` (from 27a-A8).
  - Use a dedicated `-w /nfs/scratch01/rg/fzanarello/fomo_refactor/work_g5` and outdir
    `…/fomo_refactor/g5/`.
  - Launch the head job as in the ablation study: `sbatch`, QOS `pipelines`, 1 cpu, 4 GB,
    12 h limit.
- [x] **B2. Snapshot the outputs and `compare.py` them with `tests/golden/realdata_B0.tsv`.**
  **Zero differences** is required. A difference is a stop condition (master §7c): investigate
  it, never accept it silently. (updated by D27g-3: gffcompare's order of tied entries inside a
  `.loci` line depends on the run directory, so this run may differ from B0 there alone; show it
  by re-running gffcompare on one pair in both layouts before comparing `.loci` order-insensitively.)
- [x] **B3. Record CPU-h used** in master §10. The total must stay within the budget of stop
  condition (d).

### C. Measure the goals

- [x] **C1. Fill the master §2 table with "after" values:**
  - `dev/code_metrics.py`;
  - grep counts;
  - `ruff --statistics`;
  - `nextflow lint`;
  - task counts from the V2 trace (27a-A9 against now);
  - the line counts per layer;
  - the length of CLAUDE.md.

  Every miss gets a sentence.

### D. Plan comparison and deviation validation

- [x] **D1. Diff the original against the working plan:** (updated by D27k-1) `diff -ru plans/27_refactor_original/
  plans/ --include` the `27*.md` files only. Save the raw diff to
  `/nfs/scratch01/rg/fzanarello/fomo_refactor/plan_diff.txt`.
- [x] **D2. Classify every hunk:**
  - **bookkeeping:** status, checkboxes, Results, the status log;
  - **documented deviation:** it has a `D27x-N` entry;
  - **undocumented change.**

  Every undocumented change gets a retroactive `D27x-N` entry, or is reverted if it was a
  mistake.
- [x] **D3. Validate each deviation in `plans/27z_deviation_review.md`**, in a table:

  ```
  | id | subplan | summary | class (correction / improvement / scope cut / scope addition) | consistent with U1–U8? | invariants and gate preserved? | evidence present? | verdict |
  ```

  - A deviation is **valid** only if it keeps every U-decision, keeps the gate and the
    invariants, and its evidence is reproducible (command, output or test).
  - **Invalid deviations are fixed:** change the code back toward the plan or complete the
    missing work, re-run the gates, and record it.
  - A deviation that needs a U-decision change has already been escalated (master §7). Confirm
    the user's answer is recorded.
- [x] **D4. Plan coverage:** every step of every original subplan is done, or appears in a
  validated deviation as replaced or cut. List any that were silently skipped (there must be
  none).

### E. Close-out

- [x] **E1. Write master §11 "Outcome":**
  - the goals table (before and after);
  - G5 result;
  - task-count change;
  - lines removed;
  - the deviation summary (counts by class, all valid);
  - what CI could not verify;
  - follow-ups. These include the scientific recommendations from
    [`ablation/ABLATION_REPORT.md`](https://github.com/apollo994/fomo/blob/dda0535bb62f10be53291d9a320efa46af20234c/ablation/ABLATION_REPORT.md) (adaptive donor set, per-track thresholds, mRNA coding mask,
    `splice:hq` option, decoy FDR estimator) as candidate plans 28+.
- [x] **E2. Final commit:** `refactor(27k): final validation and plan comparison`. Set the master
  status to `done (YYYY-MM-DD)`. Do not push; tell the user the branch is ready for review.

## Acceptance

- Every gate passes on the clean checkout.
- G5 has zero differences.
- `27z_deviation_review.md` has no invalid deviation.
- Master §11 is written.

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|
| D27k-1 | D1 | `diff -ru plans/27_refactor_original/ plans/ --include` the `27*.md` files | `diff -u` per file, concatenated into `plan_diff.txt` | GNU diff has no `--include` ("unrecognized option '--include=27*.md'") | none |
| D27k-2 | A1, A2 | gates on a clean checkout of `refactor` | the clean worktree is `67ccc93` (detached: `refactor` is checked out in the repository) plus the 27k fixes applied with `git apply`; `diff -r` against the working tree shows no difference but untracked tool state | D4 found work to finish (27z), and the gates must cover it; the 27k commit holds exactly that tree. G5 ran on `67ccc93` itself: of the fixes only `run_summary_tables.py` changes a code path, and it was rerun on G5's records (13/13 tables byte-identical) | none |
| D27k-3 | after E2 (CI on GitHub) | §5: other gffcompare files are compared as bytes | `.loci` files are compared with the comma-separated entries of each column sorted (`snapshot.py` rule `sorted_entries`); the golden `.loci` rows of V1–V4 and `realdata_B0.tsv` were recomputed from the original pre-refactor outputs (`fomo_refactor/golden/V*/run{1,2}`, B0) with the new rule; only `.loci` rows changed | GitHub's golden V2 failed on two mRNA `.loci` files. The same gffcompare binary (both images: sha256 `d620f7a1…`) on the same inputs writes the golden order on NFS (65,536-byte blocks) and swaps two tied query entries of one locus on local disk or tmpfs (4,096-byte blocks); CI's files are byte-identical to a local-disk rerun; `.tracking` and `.stats` never change. User decision (2026-10-16): normalise in the gate, published files unchanged | none: U1 holds (no published byte changes); master §5 has the new row |

## Results

**A. Clean-checkout gates** (worktree `fomo_refactor/clean`, env created from
`dev/environment.yml` into `fomo_refactor/env_clean`: Python 3.11.17, pytest 9.1.1, ruff 0.11.13,
pre-commit 4.7.0, gffcompare 0.12.6, gffread 0.12.7, samtools 1.21, nf-test 0.9.5, OpenJDK 21):
- G1 `pytest -q`: 277 passed.
- G2 `pre-commit run --all-files`: every hook passes (ruff, ruff format, shellcheck, nextflow lint,
  comment budget, documentation links, trailing whitespace, end of file).
- G3 `nf-test test tests/pipeline --profile +singularity`: S1–S4 pass (78 s).
- G4 `tests/golden/gate.sh fomo_refactor/gate/27k`, fresh work directories: V1, V2, V3, V4
  identical.
- Every CI `run:` command was run locally, with Singularity in place of Docker: the four above,
  plus `run_variant.sh` + `compare.py` per variant (inside `gate.sh`). The CI YAML parses.
- **Not verifiable without pushing:** the GitHub runners; `mamba-org/setup-micromamba` and its
  environment cache; `nf-core/setup-nextflow`; the `actions/cache` step for `~/.nextflow`;
  `actions/upload-artifact`; Docker execution of the pinned images, including the docker
  profile's `linux/amd64` pin for the `gff-feature-stats` modules on arm64 hosts.

**B. G5, real data.** `fomo_refactor/g5_launch/run_g5.sh` (head job 30233044, QOS `pipelines`,
1 CPU, 4 GB): `-profile crg`, the 23-species samplesheet, `--include_decoy`, fresh `-w
fomo_refactor/work_g5`, on `67ccc93`. 878 tasks succeeded in 21 min 49 s (Nextflow: 6.6 CPU-h).
- `tests/golden/compare.py tests/golden/realdata_B0.tsv fomo_refactor/g5_snapshot.tsv`:
  **7,666 expected, 7,666 actual: 0 missing, 0 extra, 0 changed.** The `.loci` files are identical
  too, so the tie-order difference D27g-3 anticipated did not occur and no order-insensitive
  comparison was needed.
- Every target has its MultiQC report (23) and the run summary is rendered.
- Task counts, B0 (pre-refactor) → G5: 1,154 → 878. Gone: GFF_STATS_TO_MQC 69,
  GFF_STATS_TARGET 46 and GFF_STATS_TARGET_TO_MQC 46 (all 23 species are `both`, U5),
  SAMTOOLS_TO_MQC 46, SEQKIT_TO_MQC 46, GFFCOMPARE_TO_SCATTER 23; GFF_STATS_PROJECTED_BATCH 46 →
  GFF_STATS_PROJECTED 46. Every other process has the same count.
- The run summary edited in 27k (D4), rerun on G5's staged records: 13/13 `summary/tables`
  files byte-identical with G5's.
- **B3.** 27k allocated 8.0 CPU-h (G5 tasks 7.2, its head job 0.4, gates and nf-test 0.4);
  cumulative 15.6 of the 40 CPU-h budget.

**C. Goals.** Measured on `67ccc93` against `c86521c` (the pre-refactor commit), in master §2.

**D. Plan comparison.** [27z_deviation_review.md](27z_deviation_review.md): 65 hunks classified;
four propagation misses fixed in the plan text; 38 deviations (two recorded retroactively, D27e-3
and D27e-4), all valid; the coverage audit found eight unfinished items in 27a, 27e, 27f, 27h and
27i, finished here, each noted in its subplan's Results.

**E.** Master §11 written; `dev/plans_index.py` rerun.

**After the push (2026-10-16), CI on GitHub.** Lint, pytest, the nf-test stub matrix (Docker) and
golden V1, V3 and V4 passed at once. Golden V2 failed on two `.loci` files
(`Lycaena_alciphron_282377.from_{Lycaena_thersamon_265360,allModels_raw}.mRNA.gffcompare.loci`):
one line of 804 in each, locus `XLOC_000201`, two query transcripts in swapped order.
- **Not the build or the image:** the gffcompare binary is identical in the Galaxy-depot
  Singularity image and the quay.io Docker image.
- **Not the directory:** rerunning the task from directories of 65–100 characters on NFS
  always gives the golden bytes.
- **The filesystem:** the same run on node-local disk or `/dev/shm` (4,096-byte blocks, NFS has
  65,536) gives CI's bytes exactly.
- `ci.yml` now uploads the snapshot and the deliverables of a failing golden variant, which is
  how CI's files were obtained.

D27k-3 follows from this. After it:
- V1–V4 of the 27k cluster gate are identical to the new golden files;
- GitHub's V2 results are identical;
- G5's 1,058 `.loci` files are identical.

On `776b1f9` every CI job passes, in both the push and
the pull-request run (GitHub Actions runs 38891360192 and 38891365163).
