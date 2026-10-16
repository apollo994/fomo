# 27a — Safety net: golden baseline, pytest, nf-test, lint, CI

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Review items: T1, T2, T3, R4.

## Status: not started

## Goal

Build the checks that every later subplan relies on **before any pipeline code changes**:

- a deterministic golden snapshot of today's deliverables;
- pytest coverage of the scripts that hold the contracts;
- an nf-test stub matrix that asserts task counts and published paths;
- lint, pre-commit and a CI workflow;
- `nextflow lint` at 0 errors.

## Scope

- **In:** `tests/`, `dev/`, `pyproject.toml`, `.pre-commit-config.yaml`, `.editorconfig`,
  `.github/workflows/ci.yml`, `nf-test.config`.
- **Code changes:** stub fixes, strict-syntax fixes and shellcheck fixes, all behaviour-neutral
  and proven by G4.
- **Out:** any other code change.

## Steps

### A. Golden baseline — on the untouched code (first commit of `refactor`)

- [ ] **A1. Write `tests/golden/snapshot.py OUTDIR > snapshot.tsv`.** It lists every
  deliverable path relative to OUTDIR with its normalised sha256, using the normalisation table
  in master §5 and these rules:
  - BAM content comes from `samtools view` (the `SAMTOOLS` env var overrides the binary).
  - Paths are sorted.
  - It also writes `snapshot.info.tsv`, which is **informational only**: for each
    `multiqc_report_data/*.txt` table, its row count and column names (used by 27i).
- [ ] **A2. Write `tests/golden/compare.py A.tsv B.tsv`.** It exits 1 on any difference and
  prints missing, extra and changed paths.
- [ ] **A3. Write `tests/golden/run_variant.sh VARIANT OUTDIR [PROFILE]`.**
  - Variants are V1–V4 (master §5). The profile defaults to `test,singularity`.
  - For V4 it builds a temporary samplesheet whose FASTAs are recompressed with plain `gzip`
    (`zcat | gzip -n`) and whose GFF3s are decompressed, under `OUTDIR/inputs/`.
  - It uses a fresh `-w OUTDIR/work`.
  - Inside SLURM it adds a `-c` file with `process.resourceLimits = [cpus: 4, memory: 16.GB,
    time: 2.h]`.
- [ ] **A4. Run V1–V4 twice each on the untouched HEAD** (T10: one `sbatch` job per run, 4 cpu,
  16 GB, ≤ 2 h). Scratch: `/nfs/scratch01/rg/fzanarello/fomo_refactor/golden/<variant>/run{1,2}`.
- [ ] **A5. Determinism check:** run 1 and run 2 of each variant must snapshot identically. For
  any file that differs:
  - find the cause;
  - only if it is real non-determinism (ordering of ties, timestamps), extend the normalisation
    and document the rule in Results;
  - otherwise stop (master §7).
- [ ] **A6. V4 must equal V1.** If it does not, record the difference as a bug for 27b. It does
  not block the baseline.
- [ ] **A7. Commit `tests/golden/V1.tsv … V4.tsv`** (from run 1), plus the informational
  `*.info.tsv` files.
- [ ] **A8. Real-data baseline.** Snapshot
  `/nfs/scratch01/rg/fzanarello/fomo_ablation/runs/B0_baseline/results` into
  `tests/golden/realdata_B0.tsv`. That directory is on scratch and is not backed up, so do it
  now. Also record the samplesheet (`/nfs/scratch01/rg/fzanarello/fomo_ablation/samplesheet_both23.csv`)
  and the command, from `runs/B0_baseline/.nextflow.log`, in `tests/golden/README.md`.
- [ ] **A9. Record the V2 trace task counts per process** (`pipeline_info/execution_trace*.txt`)
  in Results. They are the baseline for the "fewer tasks" goal and for the stub matrix.

### B. Dev environment and pytest

- [ ] **B1. Write `dev/environment.yml`** (conda-forge + bioconda) with:
  python 3.11, pytest, ruff, pre-commit, shellcheck, gffcompare 0.12.6, gffread 0.12.7,
  samtools 1.23.1 (or the nearest available; record it), and nf-test (pinned; or install the
  nf-test launcher into `dev/bin` with its version recorded). Create the env as `fomo-dev`.
- [ ] **B2. `pyproject.toml`:**
  - `[tool.pytest.ini_options]` with `testpaths = ["tests"]` and the markers `integration`
    (needs gffcompare, gffread or samtools) and `slow`;
  - `[tool.ruff]` with line length 120, select `E,F,W,B,C90,N,UP,SIM,PLR0912,PLR0915`,
    `mccabe.max-complexity = 15`, and **explicit per-file-ignores for today's offenders only**.
    Each later subplan removes entries; 27f leaves none.
- [ ] **B3. Write `tests/conftest.py`:** `BIN` path fixture, tool-availability skips, and a
  `run_script(name, *args)` helper.
- [ ] **B4. Convert `tests/curate_models/test_curate_models.py` to pytest.**
  - Keep every existing check.
  - Group them into focused `test_*` functions; shared fixtures build the model GFF3s.
  - Delete the custom `check()`/`FAILS` harness.
  - Clean up temp dirs (`tmp_path`).
- [ ] **B5. Fixture tests (new), one file per script under `tests/<script>/`:**
  - [ ] `fomo_stats`: `parse_identity` (grammar incl. `allModels_raw`, decoy), `parse_accuracy`
    on a captured `.stats` file, `f1`.
  - [ ] `gff_by_source`:
    - a split into per-source files (names from `--name-template`);
    - the duplicate transcript ID → exit 1 naming both sources;
    - empty input → no files, exit 0;
    - children inherit the source.
  - [ ] `filter_gff_by_id`:
    - keep mode with an id map + FASTA (exact match; exit 1 on a twice-mapped header, a
      missing record, or a non-transcript id);
    - drop mode with an empty or absent drop list.
  - [ ] `filter_transcript.sh` — this is the oracle for 27f, so be thorough:
    - stranded filtering (strand `?` on the transcript or a child; mixed child strands);
    - multi-exon filtering;
    - `--longest`: the spliced length decides, ties go to the smallest id, and a gene whose
      longest isoform is mono-exonic keeps its longest spliced one;
    - `Parent=a,b`;
    - a transcript without a Parent;
    - funnel counts in `--counts`;
    - plain and `.gz` input through the module's decompression path (test the script on plain
      input);
    - empty input.
  - [ ] `bam_to_gff.sh`: build a BAM from a tiny SAM with samtools. Check:
    - CIGAR with M, D, N, I and S → exon blocks;
    - single-exon models dropped by default and kept with `--include-single-exon`;
    - the `<tid>|<source>` IDs and `source=`;
    - a tid with `;` → exit 1.
  - [ ] `td2_coding_filter`:
    - `parse_coding` (complete ORF and the PSAURON threshold, `.p<N>` stripping);
    - the input stage (bare tid) against the projected stage (verbatim id);
    - an empty pep file.
  - [ ] `relocate_loci`: the same seed gives byte-identical output; decoys stay in intergenic
    intervals, keep the exon and intron shape, and are sorted by seqid; the cap is honoured.
  - [ ] `run_summary_tables`: **regression fixture.**
    - Copy the staged input directory of V1's and V2's `RUN_SUMMARY_TABLES` task, i.e. its
      `mqc/` files and roles CSV (find it with the trace hash), into
      `tests/run_summary_tables/fixtures/{V1,V2}/`. Keep it small; prune unrelated files, and
      if it exceeds about 5 MB keep V1 only.
    - Assert the outputs are byte-identical to the golden `summary/tables`.
    - 27h rewrites the input layer against this test.
  - [ ] `gff-feature-stats` smoke: `-V` runs and prints a version (skip if not x86-64 Linux).

### C. Stub run and nf-test

- [ ] **C1. Fix the stub defects** found in the review (T2):
  - [ ] the `CURATE_MODELS` stub uses `gzip`, which its image lacks. Write empty `.gz` files with
    `python3 -c "import gzip,sys; [gzip.open(p,'wt').close() for p in sys.argv[1:]]" …`;
  - [ ] the `GFF_BY_SOURCE` stub must honour the name template. Write
    `<target>.from_stub.<track>.projected.gff3`, derived from the same `ext.args` string;
  - [ ] note in `nf-test.config` that stub tests need a container profile, because the version
    `eval`s run in stub mode too.
- [ ] **C2. Install nf-test** (B1) and write `nf-test.config` (profiles `test`; `testsDir
  tests/pipeline`).
- [ ] **C3. Write `tests/pipeline/main.nf.test`** with three stub scenarios. Each asserts
  workflow success, **the task count per process**, and **the set of published paths**:
  - S1: the V1 samplesheet (F=1, D=1, S=2, T=2, T_g=1).
  - S2: `--include_mrna --include_decoy` (F=2, D=2).
  - S3: the noref samplesheet (curation without a reference, no merged file).

  Derive the expected counts from the formulas in CLAUDE.md and **cross-check them against the
  real V-run traces (A9)**. Where formula and trace disagree, the trace wins; record why.
- [ ] **C4. Pin the stub matrix as a contract.** Later subplans that change task counts on
  purpose (27i) update the expected numbers and log a deviation note in Results.

### D. Strict syntax and lint

- [ ] **D1. `nextflow lint` → 0 errors.**
  - Replace `p instanceof List ? p : [p]` (`reporting.nf:47,51`, `run_summary.nf:34`,
    `projection.nf:199`) with a v2-compatible form, e.g. `[p].flatten()`. Check that Path
    elements are not iterated: a Path is not a List, so `flatten` keeps it whole. Verify with
    `NXF_SYNTAX_PARSER=v2 nextflow run … -preview`.
  - Replace `Channel.` with `channel.`; make the implicit `it` explicit; prefix unused closure
    params with `_`.
- [ ] **D2. Make `bin/*.sh` shellcheck-clean,** behaviour-neutral (quoting and `set -euo
  pipefail` only). Leave `LC_ALL` decisions to 27b.
- [ ] **D3. Write `.editorconfig`** (spaces; 4 for nf, py and sh; LF) and
  **`.pre-commit-config.yaml`**:
  - ruff (lint only; formatting is enabled in 27f);
  - shellcheck;
  - a local hook for `nextflow lint` (errors fail; warnings are allowed until 27e);
  - trailing whitespace and end-of-file fixers on non-data files.
- [ ] **D4. Run `pre-commit run --all-files`; it must pass** with the explicit ignores from B2.

### E. CI

- [ ] **E1. Write `.github/workflows/ci.yml`** with these jobs:
  - `lint`: pre-commit;
  - `pytest`: micromamba from `dev/environment.yml`, then `pytest -q`;
  - `nf-test-stub`: docker, `nf-test test --profile +docker`;
  - `golden`: docker; V1–V4 through `run_variant.sh` with profile `test,docker`, then
    `compare.py` against `tests/golden/`.

  Pin Nextflow to 25.10.4, matching the cluster. Cache conda and the Nextflow assets.
- [ ] **E2. Validate locally:**
  - the YAML parses;
  - every `run:` command has been executed locally with Singularity in place of Docker;
  - unverifiable parts (GitHub runners, Docker) are listed in Results.

## Acceptance

- G4 baseline committed and deterministic (A5). V4 = V1, or the difference is recorded for 27b.
- G1: `pytest -q` passes.
- G2: `pre-commit run --all-files` passes.
- `nextflow lint` reports 0 errors.
- G3: `nf-test test` passes S1–S3.
- **G4 rerun after the C1/D1/D2 code changes:** V1–V4 identical to the baseline.
- CI YAML validated (E2).

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|

## Results

*(commands, numbers, determinism notes, V2 task counts, and what CI could not verify)*
