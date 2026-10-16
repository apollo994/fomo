# 27a — Safety net: golden baseline, pytest, nf-test, lint, CI

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Review items: T1, T2, T3, R4.

## Status: done (2026-10-15)

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

- [x] **A1. Write `tests/golden/snapshot.py OUTDIR > snapshot.tsv`.** It lists every
  deliverable path relative to OUTDIR with its normalised sha256, using the normalisation table
  in master §5 and these rules:
  - BAM content comes from `samtools view` (the `SAMTOOLS` env var overrides the binary).
  - Paths are sorted.
  - It also writes `snapshot.info.tsv`, which is **informational only**: for each
    `multiqc_report_data/*.txt` table, its row count and column names (used by 27i).
- [x] **A2. Write `tests/golden/compare.py A.tsv B.tsv`.** It exits 1 on any difference and
  prints missing, extra and changed paths.
- [x] **A3. Write `tests/golden/run_variant.sh VARIANT OUTDIR [PROFILE]`.**
  - Variants are V1–V4 (master §5). The profile defaults to `test,singularity`.
  - For V4 it builds a temporary samplesheet whose FASTAs are recompressed with plain `gzip`
    (`zcat | gzip -n`) and whose GFF3s are decompressed, under `OUTDIR/inputs/`.
  - It uses a fresh `-w OUTDIR/work`.
  - Inside SLURM it adds a `-c` file with `process.resourceLimits = [cpus: 4, memory: 16.GB,
    time: 2.h]`.
- [x] **A4. Run V1–V4 twice each on the untouched HEAD** (T10: one `sbatch` job per run, 4 cpu,
  16 GB, ≤ 2 h). Scratch: `/nfs/scratch01/rg/fzanarello/fomo_refactor/golden/<variant>/run{1,2}`.
- [x] **A5. Determinism check:** run 1 and run 2 of each variant must snapshot identically. For
  any file that differs:
  - find the cause;
  - only if it is real non-determinism (ordering of ties, timestamps), extend the normalisation
    and document the rule in Results;
  - otherwise stop (master §7).
- [x] **A6. V4 must equal V1.** If it does not, record the difference as a bug for 27b. It does
  not block the baseline.
- [x] **A7. Commit `tests/golden/V1.tsv … V4.tsv`** (from run 1), plus the informational
  `*.info.tsv` files.
- [x] **A8. Real-data baseline.** Snapshot
  `/nfs/scratch01/rg/fzanarello/fomo_ablation/runs/B0_baseline/results` into
  `tests/golden/realdata_B0.tsv`. That directory is on scratch and is not backed up, so do it
  now. Also record the samplesheet (`/nfs/scratch01/rg/fzanarello/fomo_ablation/samplesheet_both23.csv`)
  and the command, from `runs/B0_baseline/.nextflow.log`, in `tests/golden/README.md`.
- [x] **A9. Record the V2 trace task counts per process** (`pipeline_info/execution_trace*.txt`)
  in Results. They are the baseline for the "fewer tasks" goal and for the stub matrix.

### B. Dev environment and pytest

- [x] **B1. Write `dev/environment.yml`** (conda-forge + bioconda) with:
  python 3.11, pytest, ruff, pre-commit, shellcheck, gffcompare 0.12.6, gffread 0.12.7,
  samtools 1.23.1 (or the nearest available; record it), and nf-test (pinned; or install the
  nf-test launcher into `dev/bin` with its version recorded). Create the env as `fomo-dev`.
- [x] **B2. `pyproject.toml`:**
  - `[tool.pytest.ini_options]` with `testpaths = ["tests"]` and the markers `integration`
    (needs gffcompare, gffread or samtools) and `slow`;
  - `[tool.ruff]` with line length 120, select `E,F,W,B,C90,N,UP,SIM,PLR0912,PLR0915`,
    `mccabe.max-complexity = 15`, and **explicit per-file-ignores for today's offenders only**.
    Each later subplan removes entries; 27f leaves none.
- [x] **B3. Write `tests/conftest.py`:** `BIN` path fixture, tool-availability skips, and a
  `run_script(name, *args)` helper.
- [x] **B4. Convert `tests/curate_models/test_curate_models.py` to pytest.**
  - Keep every existing check.
  - Group them into focused `test_*` functions; shared fixtures build the model GFF3s.
  - Delete the custom `check()`/`FAILS` harness.
  - Clean up temp dirs (`tmp_path`).
- [x] **B5. Fixture tests (new), one file per script under `tests/<script>/`:**
  - [x] `fomo_stats`: `parse_identity` (grammar incl. `allModels_raw`, decoy), `parse_accuracy`
    on a captured `.stats` file, `f1`.
  - [x] `gff_by_source`:
    - a split into per-source files (names from `--name-template`);
    - the duplicate transcript ID → exit 1 naming both sources;
    - empty input → no files, exit 0;
    - children inherit the source.
  - [x] `filter_gff_by_id`:
    - keep mode with an id map + FASTA (exact match; exit 1 on a twice-mapped header, a
      missing record, or a non-transcript id);
    - drop mode with an empty or absent drop list.
  - [x] `filter_transcript.sh` — this is the oracle for 27f, so be thorough:
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
  - [x] `bam_to_gff.sh`: build a BAM from a tiny SAM with samtools. Check:
    - CIGAR with M, D, N, I and S → exon blocks;
    - single-exon models dropped by default and kept with `--include-single-exon`;
    - the `<tid>|<source>` IDs and `source=`;
    - a tid with `;` → exit 1.
  - [x] `td2_coding_filter`:
    - `parse_coding` (complete ORF and the PSAURON threshold, `.p<N>` stripping);
    - the input stage (bare tid) against the projected stage (verbatim id);
    - an empty pep file.
  - [x] `relocate_loci`: the same seed gives byte-identical output; decoys stay in intergenic
    intervals, keep the exon and intron shape, and are sorted by seqid; the cap is honoured.
  - [x] `run_summary_tables`: **regression fixture.**
    - Copy the staged input directory of V1's and V2's `RUN_SUMMARY_TABLES` task, i.e. its
      `mqc/` files and roles CSV (find it with the trace hash), into
      `tests/run_summary_tables/fixtures/{V1,V2}/`. Keep it small; prune unrelated files, and
      if it exceeds about 5 MB keep V1 only.
    - Assert the outputs are byte-identical to the golden `summary/tables`.
    - 27h rewrites the input layer against this test.
  - [x] `gff-feature-stats` smoke: `-V` runs and prints a version (skip if not x86-64 Linux).

### C. Stub run and nf-test

- [x] **C1. Fix the stub defects** found in the review (T2):
  - [x] the `CURATE_MODELS` stub uses `gzip`, which its image lacks. Write empty `.gz` files with
    `python3 -c "import gzip,sys; [gzip.open(p,'wt').close() for p in sys.argv[1:]]" …`;
  - [x] the `GFF_BY_SOURCE` stub must honour the name template. Write
    `<target>.from_stub.<track>.projected.gff3`, derived from the same `ext.args` string;
  - [x] note in `nf-test.config` that stub tests need a container profile, because the version
    `eval`s run in stub mode too.
- [x] **C2. Install nf-test** (B1) and write `nf-test.config` (profiles `test`; `testsDir
  tests/pipeline`).
- [x] **C3. Write `tests/pipeline/main.nf.test`** with three stub scenarios. Each asserts
  workflow success, **the task count per process**, and **the set of published paths**:
  - S1: the V1 samplesheet (F=1, D=1, S=2, T=2, T_g=1).
  - S2: `--include_mrna --include_decoy` (F=2, D=2).
  - S3: the noref samplesheet (curation without a reference, no merged file).

  Derive the expected counts from the formulas in CLAUDE.md and **cross-check them against the
  real V-run traces (A9)**. Where formula and trace disagree, the trace wins; record why.
- [x] **C4. Pin the stub matrix as a contract.** Later subplans that change task counts on
  purpose (27i) update the expected numbers and log a deviation note in Results.

### D. Strict syntax and lint

- [x] **D1. `nextflow lint` → 0 errors.**
  - Replace `p instanceof List ? p : [p]` (`reporting.nf:47,51`, `run_summary.nf:34`,
    `projection.nf:199`) with a v2-compatible form, e.g. `[p].flatten()`. Check that Path
    elements are not iterated: a Path is not a List, so `flatten` keeps it whole. Verify with
    `NXF_SYNTAX_PARSER=v2 nextflow run … -preview`.
  - Replace `Channel.` with `channel.`; make the implicit `it` explicit; prefix unused closure
    params with `_`.
- [x] **D2. Make `bin/*.sh` shellcheck-clean,** behaviour-neutral (quoting and `set -euo
  pipefail` only). Leave `LC_ALL` decisions to 27b.
- [x] **D3. Write `.editorconfig`** (spaces; 4 for nf, py and sh; LF) and
  **`.pre-commit-config.yaml`**:
  - ruff (lint only; formatting is enabled in 27f);
  - shellcheck;
  - a local hook for `nextflow lint` (errors fail; warnings are allowed until 27e);
  - trailing whitespace and end-of-file fixers on non-data files.
- [x] **D4. Run `pre-commit run --all-files`; it must pass** with the explicit ignores from B2.

### E. CI

- [x] **E1. Write `.github/workflows/ci.yml`** with these jobs:
  - `lint`: pre-commit;
  - `pytest`: micromamba from `dev/environment.yml`, then `pytest -q`;
  - `nf-test-stub`: docker, `nf-test test --profile +docker`;
  - `golden`: docker; V1–V4 through `run_variant.sh` with profile `test,docker`, then
    `compare.py` against `tests/golden/`.

  Pin Nextflow to 25.10.4, matching the cluster. Cache conda and the Nextflow assets.
- [x] **E2. Validate locally:**
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
| D27a-1 | A (tooling) | G4 is run by hand: `sbatch` each `run_variant.sh`, then `compare.py` | added `tests/golden/gate.sh SCRATCH`, which submits V1–V4 as SLURM jobs, waits, and compares all four | every later subplan runs G4; one command removes four manual submissions and comparisons | later subplans run G4 as `tests/golden/gate.sh <scratch>/gate/27x` |
| D27a-2 | B2, D3 | ruff formatting is enabled in 27f | `ruff format` (quote style preserved) is configured now and applied **to the new test files only**; `bin/` is formatted in 27f as planned. Lint debt also covers UP007/UP038/B905, which ruff reports once `target-version = "py311"` is set | new files start formatted; the extra codes are today's offenders, not new debt | 27f removes these ignores together with the others |
| D27a-3 | C3 | the stub matrix asserts the published paths | the expected path lists come from the stub runs. They are a subset of the real runs' paths: the run-summary stub writes two tables, the gffcompare stub no `.annotated.gtf`/`.tmap`/`.refmap`, and per-source names read `from_stub` | stubs cannot know source names; the real path set is guarded by G4 | when a subplan changes stubs or published names, regenerate `tests/pipeline/expected/*.paths` from a stub run and review the diff |
| D27a-4 | (finding) | — | the V1 run logs `There's no process matching config selector: …EXTRACT_DECOY_SEQUENCES` (and the other decoy processes) when `--include_decoy` is off | evidence that `if`-gating an *included* process does raise the selector warning | 27e-C2 is already answered: keep `ext.when` for curation, or make the selectors anchored `.*:NAME$` for gated processes, and re-check the warning (updated 27e-C2) |

## Results

> **Correction (found in 27k D4):** D3's end-of-file fixer and E1's Nextflow-assets cache were
> not made, and B1 left nf-test unpinned; Results below did not say so. 27k added
> `dev/end_of_file.py` as a pre-commit hook (tested in `tests/dev`), an `actions/cache` step for
> `~/.nextflow` in the three CI jobs that run Nextflow, and `nf-test=0.9.5` in `dev/environment.yml`.

**Golden baseline (A).**
- Untouched code at `c86521c`. Each variant was run twice as a SLURM job (4 cpu, test profile,
  singularity). Every run took about 1–2 min.
- **Deterministic:** run 1 and run 2 are identical for V1–V4, with no normalisation beyond §5.
- **V4 = V1:** plain-gzip FASTA and uncompressed GFF3 give identical deliverables.
- Deliverables per variant: V1 36, V2 119, V3 35, V4 36. Stored as `tests/golden/V*.tsv`, plus
  `V*.info.tsv` (61 MultiQC data tables each, informational).
- **Real data:** `realdata_B0.tsv` has 7,666 deliverables (36 BAMs, 805 gffcompare stats, …)
  of the 23-species lncRNA+decoy run. The samplesheet and command are in `tests/golden/README.md`.

**Task counts of the real runs (A9).** These are the baseline for the "fewer tasks" goal. Totals
are V1 = 61, V2 = 167, V3 = 59.

V2 per process:
- **Batch and filter steps:** BAM_TO_GFF 8, CURATE_MODELS 8, FILTER_ALLMODELS 8,
  FILTER_ANNOTATION 4, FILTER_TRANSCRIPT 4, MINIMAP2_ALIGN 8, MINIMAP2_INDEX 2,
  MERGE_SOURCE_FASTA 4, SPLIT_GFF_BY_SOURCE 8, RENAME_FASTA_HEADERS 8.
- **Decoy and extraction steps:** BEDTOOLS_COMPLEMENT 2, FASTA_SIZES 2, GFF_TO_GENE_BED 2,
  RELOCATE_LOCI 4, EXTRACT_DECOY_SEQUENCES 4, EXTRACT_PROJECTED_LNC 2, EXTRACT_SEQUENCES 4.
- **Stats and adapters:** GFF_STATS 10, GFF_STATS_TO_MQC 10, GFF_STATS_TARGET 4,
  GFF_STATS_TARGET_TO_MQC 4, GFF_STATS_PROJECTED_BATCH 8, SAMTOOLS_STATS 8, SAMTOOLS_TO_MQC 8,
  SEQKIT_STATS 8, SEQKIT_TO_MQC 8, TD2_PREDICT 4, TD2_CODING_FILTER 4.
- **Benchmarking and reports:** GFFCOMPARE_BATCH 4, GFFCOMPARE_TO_SCATTER 1, MULTIQC 2,
  MULTIQC_RUN_SUMMARY 1, RUN_SUMMARY_TABLES 1.

The stub runs reproduce exactly these counts.

**Tests (B).**
- **Dev env:** `dev/environment.yml` → `fomo-dev` with python 3.11, pytest 9, ruff 0.11,
  shellcheck 0.11, gffcompare 0.12.6, gffread 0.12.7, nf-test 0.9.5 and openjdk 21.
  - samtools is 1.21. The pipeline's image uses 1.23.1, but the dev env only needs samtools to
    build test BAMs and read BAM records.
- **pytest: 81 tests pass.**
  - The curate suite was converted: 23 tests, every original check kept.
  - New: fomo_stats 13, gff_by_source 6, filter_gff_by_id 7, filter_transcript 14 (the oracle
    for 27f), bam_to_gff 6, td2_coding_filter 4, relocate_loci 4, gff-feature-stats 2.
  - **run_summary_tables 2:** byte-identical regression against V1/V2's published tables, from
    the staged inputs pruned to the 9 suffixes the script reads.

**Stub matrix (C).**
- Stub defects fixed: the CURATE_MODELS stub no longer calls `gzip` (absent from its image), and
  the GFF_BY_SOURCE stub honours `--name-template`.
- `nf-test test tests/pipeline --profile +singularity`: S1–S3 pass (task counts per process +
  published paths; about 75 s, as a SLURM job).

**Lint (D).**
- `nextflow lint`: 0 errors and **0 warnings**. Fixed `[p].flatten()` ×4, `channel.` ×5,
  explicit closure params ×3 and unused params ×7.
- `NXF_SYNTAX_PARSER=v2` handles the `flatten` idiom.
- shellcheck was already clean.
- `pre-commit run --all-files` passes: ruff, shellcheck, nextflow lint, trailing whitespace.

**CI (E).** `.github/workflows/ci.yml` (lint, pytest, nf-test-stub, golden V1–V4) parses. Every
command was run locally with Singularity in place of Docker.
- **Unverified until the branch is pushed:** the GitHub runners, the micromamba environment
  caching, Docker execution, and the `nf-core/setup-nextflow` action.

**Gate after the 27a code changes:** `tests/golden/gate.sh …/gate/27a` → V1–V4 identical.

**Cluster use:** about 1.5 CPU-h (golden ×8, gate ×4, real-data snapshot, nf-test).
